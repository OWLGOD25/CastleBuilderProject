/** @file Week4-6-ShapeComplete.cpp
 *  @brief Shape Practice Solution.
 *
 *  Place all of the scene geometry in one big vertex and index buffer.
 * Then use the DrawIndexedInstanced method to draw one object at a time ((as the
 * world matrix needs to be changed between objects)
 *
 *   Controls:
 *   Hold down '1' key to view scene in wireframe mode.
 *   Hold the left mouse button down and move the mouse to rotate.
 *   Hold the right mouse button down and move the mouse to zoom in and out.
 *
 *  @author Hooman Salamat
 */


#include "../../Common/d3dApp.h"
#include "../../Common/MathHelper.h"
#include "../../Common/UploadBuffer.h"
#include "../../Common/GeometryGenerator.h"
#include "FrameResource.h"
#include "../../Common/DDSTextureLoader.h"
#include <array>
#include <filesystem>
#include <string>


using Microsoft::WRL::ComPtr;
using namespace DirectX;
using namespace DirectX::PackedVector;

const int gNumFrameResources = 3;

struct RenderItem
{
	RenderItem() = default;

	XMFLOAT4X4 World = MathHelper::Identity4x4();
	XMFLOAT4X4 TexTransform = MathHelper::Identity4x4();

	int NumFramesDirty = gNumFrameResources;
	UINT ObjCBIndex = -1;

	Material* Mat = nullptr;
	MeshGeometry* Geo = nullptr;

	D3D12_PRIMITIVE_TOPOLOGY PrimitiveType = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;

	UINT IndexCount = 0;
	UINT StartIndexLocation = 0;
	int BaseVertexLocation = 0;
};

class ShapesApp : public D3DApp
{
public:
	ShapesApp(HINSTANCE hInstance);
	ShapesApp(const ShapesApp& rhs) = delete;
	ShapesApp& operator=(const ShapesApp& rhs) = delete;
	~ShapesApp();

	virtual bool Initialize()override;

private:
	virtual void OnResize()override;
	virtual void Update(const GameTimer& gt)override;
	virtual void Draw(const GameTimer& gt)override;

	virtual void OnMouseDown(WPARAM btnState, int x, int y)override;
	virtual void OnMouseUp(WPARAM btnState, int x, int y)override;
	virtual void OnMouseMove(WPARAM btnState, int x, int y)override;

	void OnKeyboardInput(const GameTimer& gt);
	void UpdateCamera(const GameTimer& gt);
	void AnimateMaterials(const GameTimer& gt);
	void UpdateObjectCBs(const GameTimer& gt);
	void UpdateMainPassCB(const GameTimer& gt);

	void BuildDescriptorHeaps();
	void BuildConstantBufferViews();
	void BuildRootSignature();
	void BuildShadersAndInputLayout();
	void BuildShapeGeometry();
	void BuildPSOs();
	void BuildFrameResources();
	void BuildRenderItems();
	void DrawRenderItems(ID3D12GraphicsCommandList* cmdList, const std::vector<RenderItem*>& ritems);
	// NEW: returns default samplers for textures
	std::array<const CD3DX12_STATIC_SAMPLER_DESC, 6> GetStaticSamplers();

	void LoadTextures();
	void BuildMaterials();
	void UpdateMaterialCBs(const GameTimer& gt);

private:

	std::vector<std::unique_ptr<FrameResource>> mFrameResources;
	FrameResource* mCurrFrameResource = nullptr;
	int mCurrFrameResourceIndex = 0;

	ComPtr<ID3D12RootSignature> mRootSignature = nullptr;
	ComPtr<ID3D12DescriptorHeap> mCbvHeap = nullptr;

	ComPtr<ID3D12DescriptorHeap> mSrvDescriptorHeap = nullptr;

	std::unordered_map<std::string, std::unique_ptr<MeshGeometry>> mGeometries;
	std::unordered_map<std::string, ComPtr<ID3DBlob>> mShaders;
	std::unordered_map<std::string, ComPtr<ID3D12PipelineState>> mPSOs;

	// NEW: stores all textures (stone, wood, grass)
	std::unordered_map<std::string, std::unique_ptr<Texture>> mTextures;
	// NEW: stores all materials like stone, wood, grass
	std::unordered_map<std::string, std::unique_ptr<Material>> mMaterials;

	std::vector<D3D12_INPUT_ELEMENT_DESC> mInputLayout;

	// List of all the render items.
	std::vector<std::unique_ptr<RenderItem>> mAllRitems;

	// Render items divided by PSO.
	std::vector<RenderItem*> mOpaqueRitems;

	std::vector<RenderItem*> mTransparentRitems;

	PassConstants mMainPassCB;

	UINT mPassCbvOffset = 0;

	bool mIsWireframe = false;

	XMFLOAT3 mEyePos = { 0.0f, 0.0f, 0.0f };
	XMFLOAT4X4 mView = MathHelper::Identity4x4();
	XMFLOAT4X4 mProj = MathHelper::Identity4x4();

	float mTheta = 1.5f * XM_PI;
	float mPhi = 0.2f * XM_PI;
	float mRadius = 15.0f;

	POINT mLastMousePos;
};

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE prevInstance,
	PSTR cmdLine, int showCmd)
{
	// Enable run-time memory check for debug builds.
#if defined(DEBUG) | defined(_DEBUG)
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

	try
	{
		ShapesApp theApp(hInstance);
		if (!theApp.Initialize())
			return 0;

		return theApp.Run();
	}
	catch (DxException& e)
	{
		MessageBox(nullptr, e.ToString().c_str(), L"HR Failed", MB_OK);
		return 0;
	}
}

ShapesApp::ShapesApp(HINSTANCE hInstance)
	: D3DApp(hInstance)
{
}

ShapesApp::~ShapesApp()
{
	if (md3dDevice != nullptr)
		FlushCommandQueue();
}

bool ShapesApp::Initialize()
{
	if (!D3DApp::Initialize())
		return false;

	// Reset the command list to prep for initialization commands.
	ThrowIfFailed(mCommandList->Reset(mDirectCmdListAlloc.Get(), nullptr));

	LoadTextures();              // NEW: load texture files first
	BuildRootSignature();
	BuildShadersAndInputLayout();
	BuildShapeGeometry();
	BuildMaterials();            // NEW: create stone/wood/grass materials
	BuildRenderItems();
	BuildFrameResources();
	BuildDescriptorHeaps();
	BuildConstantBufferViews();
	BuildPSOs();

	// Execute the initialization commands.
	ThrowIfFailed(mCommandList->Close());
	ID3D12CommandList* cmdsLists[] = { mCommandList.Get() };
	mCommandQueue->ExecuteCommandLists(_countof(cmdsLists), cmdsLists);

	// Wait until initialization is complete.
	FlushCommandQueue();

	return true;
}


// ============================================================
// here is to add the textrues
// ============================================================
void ShapesApp::LoadTextures()
{
	auto stone1Tex = std::make_unique<Texture>();
	stone1Tex->Name = "stone1Tex";
	stone1Tex->Filename = L"Textures/stone1.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(
		md3dDevice.Get(),
		mCommandList.Get(),
		stone1Tex->Filename.c_str(),
		stone1Tex->Resource,
		stone1Tex->UploadHeap));
	mTextures["stone1Tex"] = std::move(stone1Tex);

	auto woodTex = std::make_unique<Texture>();
	woodTex->Name = "woodTex";
	woodTex->Filename = L"Textures/wood.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(
		md3dDevice.Get(),
		mCommandList.Get(),
		woodTex->Filename.c_str(),
		woodTex->Resource,
		woodTex->UploadHeap));
	mTextures["woodTex"] = std::move(woodTex);

	auto grassTex = std::make_unique<Texture>();
	grassTex->Name = "grassTex";
	grassTex->Filename = L"Textures/grass.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(
		md3dDevice.Get(),
		mCommandList.Get(),
		grassTex->Filename.c_str(),
		grassTex->Resource,
		grassTex->UploadHeap));
	mTextures["grassTex"] = std::move(grassTex);

	auto roofTex = std::make_unique<Texture>();
	roofTex->Name = "roofTex";
	roofTex->Filename = L"Textures/roof.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(
		md3dDevice.Get(),
		mCommandList.Get(),
		roofTex->Filename.c_str(),
		roofTex->Resource,
		roofTex->UploadHeap));
	mTextures["roofTex"] = std::move(roofTex);

	auto towerTex = std::make_unique<Texture>();
	towerTex->Name = "towerTex";
	towerTex->Filename = L"Textures/tower.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(
		md3dDevice.Get(),
		mCommandList.Get(),
		towerTex->Filename.c_str(),
		towerTex->Resource,
		towerTex->UploadHeap));
	mTextures["towerTex"] = std::move(towerTex);

	auto triprismTex = std::make_unique<Texture>();
	triprismTex->Name = "triprismTex";
	triprismTex->Filename = L"Textures/triprism.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(
		md3dDevice.Get(),
		mCommandList.Get(),
		triprismTex->Filename.c_str(),
		triprismTex->Resource,
		triprismTex->UploadHeap));
	mTextures["triprismTex"] = std::move(triprismTex);

	auto torusTex = std::make_unique<Texture>();
	torusTex->Name = "torusTex";
	torusTex->Filename = L"Textures/torus.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(
		md3dDevice.Get(),
		mCommandList.Get(),
		torusTex->Filename.c_str(),
		torusTex->Resource,
		torusTex->UploadHeap));
	mTextures["torusTex"] = std::move(torusTex);

	auto diamondTex = std::make_unique<Texture>();
	diamondTex->Name = "diamondTex";
	diamondTex->Filename = L"Textures/diamond.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(
		md3dDevice.Get(),
		mCommandList.Get(),
		diamondTex->Filename.c_str(),
		diamondTex->Resource,
		diamondTex->UploadHeap));
	mTextures["diamondTex"] = std::move(diamondTex);

	auto waterTex = std::make_unique<Texture>();
	waterTex->Name = "waterTex";
	waterTex->Filename = L"Textures/water1.dds";
	ThrowIfFailed(DirectX::CreateDDSTextureFromFile12(
		md3dDevice.Get(),
		mCommandList.Get(),
		waterTex->Filename.c_str(),
		waterTex->Resource,
		waterTex->UploadHeap));
	mTextures["waterTex"] = std::move(waterTex);
}

void ShapesApp::OnResize()
{
	D3DApp::OnResize();

	// The window resized, so update the aspect ratio and recompute the projection matrix.
	XMMATRIX P = XMMatrixPerspectiveFovLH(0.25f * MathHelper::Pi, AspectRatio(), 1.0f, 1000.0f);
	XMStoreFloat4x4(&mProj, P);
}

void ShapesApp::Update(const GameTimer& gt)
{
	OnKeyboardInput(gt);
	UpdateCamera(gt);

	// Cycle through the circular frame resource array.
	mCurrFrameResourceIndex = (mCurrFrameResourceIndex + 1) % gNumFrameResources;
	mCurrFrameResource = mFrameResources[mCurrFrameResourceIndex].get();

	// Has the GPU finished processing the commands of the current frame resource?
	// If not, wait until the GPU has completed commands up to this fence point.
	if (mCurrFrameResource->Fence != 0 && mFence->GetCompletedValue() < mCurrFrameResource->Fence)
	{
		HANDLE eventHandle = CreateEventEx(nullptr, nullptr, false, EVENT_ALL_ACCESS);
		ThrowIfFailed(mFence->SetEventOnCompletion(mCurrFrameResource->Fence, eventHandle));
		WaitForSingleObject(eventHandle, INFINITE);
		CloseHandle(eventHandle);
	}

	AnimateMaterials(gt);
	UpdateObjectCBs(gt);
	UpdateMaterialCBs(gt);   // NEW: update materials too
	UpdateMainPassCB(gt);
}

// ============================================================
// PART 3: WATER RENDER ORDER
// Opaque objects are drawn first, then the water is drawn
// using the transparent PSO so blending works correctly.
// ============================================================
void ShapesApp::Draw(const GameTimer& gt)
{
	auto cmdListAlloc = mCurrFrameResource->CmdListAlloc;

	// Reuse the memory associated with command recording.
	// We can only reset when the associated command lists have finished execution on the GPU.
	ThrowIfFailed(cmdListAlloc->Reset());

	// A command list can be reset after it has been added to the command queue via ExecuteCommandList.
	// Reusing the command list reuses memory.
	if (mIsWireframe)
	{
		ThrowIfFailed(mCommandList->Reset(cmdListAlloc.Get(), mPSOs["opaque_wireframe"].Get()));
	}
	else
	{
		ThrowIfFailed(mCommandList->Reset(cmdListAlloc.Get(), mPSOs["opaque"].Get()));
	}

	mCommandList->RSSetViewports(1, &mScreenViewport);
	mCommandList->RSSetScissorRects(1, &mScissorRect);

	// Indicate a state transition on the resource usage.
	mCommandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(CurrentBackBuffer(),
		D3D12_RESOURCE_STATE_PRESENT, D3D12_RESOURCE_STATE_RENDER_TARGET));

	// Clear the back buffer and depth buffer.
	mCommandList->ClearRenderTargetView(CurrentBackBufferView(), Colors::LightSteelBlue, 0, nullptr);
	mCommandList->ClearDepthStencilView(DepthStencilView(), D3D12_CLEAR_FLAG_DEPTH | D3D12_CLEAR_FLAG_STENCIL, 1.0f, 0, 0, nullptr);

	// Specify the buffers we are going to render to.
	mCommandList->OMSetRenderTargets(1, &CurrentBackBufferView(), true, &DepthStencilView());

	ID3D12DescriptorHeap* descriptorHeaps[] = { mSrvDescriptorHeap.Get() };
	mCommandList->SetDescriptorHeaps(_countof(descriptorHeaps), descriptorHeaps);

	mCommandList->SetGraphicsRootSignature(mRootSignature.Get());

	auto passCB = mCurrFrameResource->PassCB->Resource();
	mCommandList->SetGraphicsRootConstantBufferView(2, passCB->GetGPUVirtualAddress());

	DrawRenderItems(mCommandList.Get(), mOpaqueRitems);

	mCommandList->SetPipelineState(mPSOs["transparent"].Get());
	DrawRenderItems(mCommandList.Get(), mTransparentRitems);

	// Indicate a state transition on the resource usage.
	mCommandList->ResourceBarrier(1, &CD3DX12_RESOURCE_BARRIER::Transition(CurrentBackBuffer(),
		D3D12_RESOURCE_STATE_RENDER_TARGET, D3D12_RESOURCE_STATE_PRESENT));

	// Done recording commands.
	ThrowIfFailed(mCommandList->Close());

	// Add the command list to the queue for execution.
	ID3D12CommandList* cmdsLists[] = { mCommandList.Get() };
	mCommandQueue->ExecuteCommandLists(_countof(cmdsLists), cmdsLists);

	// Swap the back and front buffers
	ThrowIfFailed(mSwapChain->Present(0, 0));
	mCurrBackBuffer = (mCurrBackBuffer + 1) % SwapChainBufferCount;

	// Advance the fence value to mark commands up to this fence point.
	mCurrFrameResource->Fence = ++mCurrentFence;

	// Add an instruction to the command queue to set a new fence point. 
	// Because we are on the GPU timeline, the new fence point won't be 
	// set until the GPU finishes processing all the commands prior to this Signal().
	mCommandQueue->Signal(mFence.Get(), mCurrentFence);
}

void ShapesApp::OnMouseDown(WPARAM btnState, int x, int y)
{
	mLastMousePos.x = x;
	mLastMousePos.y = y;

	SetCapture(mhMainWnd);
}

void ShapesApp::OnMouseUp(WPARAM btnState, int x, int y)
{
	ReleaseCapture();
}

void ShapesApp::OnMouseMove(WPARAM btnState, int x, int y)
{
	if ((btnState & MK_LBUTTON) != 0)
	{
		// Make each pixel correspond to a quarter of a degree.
		float dx = XMConvertToRadians(0.25f * static_cast<float>(x - mLastMousePos.x));
		float dy = XMConvertToRadians(0.25f * static_cast<float>(y - mLastMousePos.y));

		// Update angles based on input to orbit camera around box.
		mTheta += dx;
		mPhi += dy;

		// Restrict the angle mPhi.
		mPhi = MathHelper::Clamp(mPhi, 0.1f, MathHelper::Pi - 0.1f);
	}
	else if ((btnState & MK_RBUTTON) != 0)
	{
		// Make each pixel correspond to 0.2 unit in the scene.
		float dx = 0.05f * static_cast<float>(x - mLastMousePos.x);
		float dy = 0.05f * static_cast<float>(y - mLastMousePos.y);

		// Update the camera radius based on input.
		mRadius += dx - dy;

		// Restrict the radius.
		mRadius = MathHelper::Clamp(mRadius, 5.0f, 150.0f);
	}

	mLastMousePos.x = x;
	mLastMousePos.y = y;
}

void ShapesApp::OnKeyboardInput(const GameTimer& gt)
{
	if (GetAsyncKeyState('1') & 0x8000)
		mIsWireframe = true;
	else
		mIsWireframe = false;
}

void ShapesApp::UpdateCamera(const GameTimer& gt)
{
	// Convert Spherical to Cartesian coordinates.
	mEyePos.x = mRadius * sinf(mPhi) * cosf(mTheta);
	mEyePos.z = mRadius * sinf(mPhi) * sinf(mTheta);
	mEyePos.y = mRadius * cosf(mPhi);

	// Build the view matrix.
	XMVECTOR pos = XMVectorSet(mEyePos.x, mEyePos.y, mEyePos.z, 1.0f);
	XMVECTOR target = XMVectorZero();
	XMVECTOR up = XMVectorSet(0.0f, 1.0f, 0.0f, 0.0f);

	XMMATRIX view = XMMatrixLookAtLH(pos, target, up);
	XMStoreFloat4x4(&mView, view);
}

// ============================================================
// PART 3: WATER ANIMATION
// Scrolls the water material texture coordinates over time
// to make the water look like it is moving.
// ============================================================
void ShapesApp::AnimateMaterials(const GameTimer& gt)
{
	auto waterMat = mMaterials["water"].get();

	// Access correctly using ._41 and ._42 (translation part)
	float& tu = waterMat->MatTransform._41;
	float& tv = waterMat->MatTransform._42;

	tu += 0.1f * gt.DeltaTime();
	tv += 0.02f * gt.DeltaTime();

	if (tu >= 1.0f) tu -= 1.0f;
	if (tv >= 1.0f) tv -= 1.0f;

	waterMat->NumFramesDirty = gNumFrameResources;
}

void ShapesApp::UpdateObjectCBs(const GameTimer& gt)
{
	auto currObjectCB = mCurrFrameResource->ObjectCB.get();

	for (auto& e : mAllRitems)
	{
		if (e->NumFramesDirty > 0)
		{
			XMMATRIX world = XMLoadFloat4x4(&e->World);
			XMMATRIX texTransform = XMLoadFloat4x4(&e->TexTransform);

			ObjectConstants objConstants;
			XMStoreFloat4x4(&objConstants.World, XMMatrixTranspose(world));
			XMStoreFloat4x4(&objConstants.TexTransform, XMMatrixTranspose(texTransform));

			currObjectCB->CopyData(e->ObjCBIndex, objConstants);
			e->NumFramesDirty--;
		}
	}
}

// ============================================================
// PART 3: MATERIAL BUFFER UPDATE
// Sends updated material data, including the animated water
// material transform, to the GPU every frame.
// ============================================================
void ShapesApp::UpdateMaterialCBs(const GameTimer& gt)
{
	auto currMaterialCB = mCurrFrameResource->MaterialCB.get();

	for (auto& e : mMaterials)
	{
		auto mat = e.second.get();

		if (mat->NumFramesDirty > 0)
		{
			XMMATRIX matTransform = XMLoadFloat4x4(&mat->MatTransform);

			MaterialConstants matConstants;
			matConstants.DiffuseAlbedo = mat->DiffuseAlbedo;
			matConstants.FresnelR0 = mat->FresnelR0;
			matConstants.Roughness = mat->Roughness;
			XMStoreFloat4x4(&matConstants.MatTransform, XMMatrixTranspose(matTransform));

			currMaterialCB->CopyData(mat->MatCBIndex, matConstants);

			mat->NumFramesDirty--;
		}
	}
}

void ShapesApp::UpdateMainPassCB(const GameTimer& gt)
{
	XMMATRIX view = XMLoadFloat4x4(&mView);
	XMMATRIX proj = XMLoadFloat4x4(&mProj);

	XMMATRIX viewProj = XMMatrixMultiply(view, proj);
	XMMATRIX invView = XMMatrixInverse(&XMMatrixDeterminant(view), view);
	XMMATRIX invProj = XMMatrixInverse(&XMMatrixDeterminant(proj), proj);
	XMMATRIX invViewProj = XMMatrixInverse(&XMMatrixDeterminant(viewProj), viewProj);

	XMStoreFloat4x4(&mMainPassCB.View, XMMatrixTranspose(view));
	XMStoreFloat4x4(&mMainPassCB.InvView, XMMatrixTranspose(invView));
	XMStoreFloat4x4(&mMainPassCB.Proj, XMMatrixTranspose(proj));
	XMStoreFloat4x4(&mMainPassCB.InvProj, XMMatrixTranspose(invProj));
	XMStoreFloat4x4(&mMainPassCB.ViewProj, XMMatrixTranspose(viewProj));
	XMStoreFloat4x4(&mMainPassCB.InvViewProj, XMMatrixTranspose(invViewProj));
	mMainPassCB.EyePosW = mEyePos;
	mMainPassCB.RenderTargetSize = XMFLOAT2((float)mClientWidth, (float)mClientHeight);
	mMainPassCB.InvRenderTargetSize = XMFLOAT2(1.0f / mClientWidth, 1.0f / mClientHeight);
	mMainPassCB.gAmbientLight = XMFLOAT4(0.25f, 0.2f, 0.15f, 1.0f);
	mMainPassCB.gLightDir = XMFLOAT3(-0.6f, -0.7f, 0.3f);
	mMainPassCB.gLightStrength = XMFLOAT3(1.0f, 0.85f, 0.6f);
	mMainPassCB.NearZ = 1.0f;
	mMainPassCB.FarZ = 1000.0f;
	mMainPassCB.TotalTime = gt.TotalTime();
	mMainPassCB.DeltaTime = gt.DeltaTime();

	auto currPassCB = mCurrFrameResource->PassCB.get();
	currPassCB->CopyData(0, mMainPassCB);
}


// ============================================================
// PART 3: WATER TEXTURE SRV SETUP
// Creates the SRV heap for all scene textures and assigns the
// water texture to its descriptor slot so it can be sampled
// by the shader.
// ============================================================
void ShapesApp::BuildDescriptorHeaps()
{
	// NEW: create SRV heap for 9 textures
	D3D12_DESCRIPTOR_HEAP_DESC srvHeapDesc = {};
	srvHeapDesc.NumDescriptors = 9;
	srvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
	srvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
	srvHeapDesc.NodeMask = 0;
	ThrowIfFailed(md3dDevice->CreateDescriptorHeap(&srvHeapDesc, IID_PPV_ARGS(&mSrvDescriptorHeap)));

	CD3DX12_CPU_DESCRIPTOR_HANDLE hDescriptor(mSrvDescriptorHeap->GetCPUDescriptorHandleForHeapStart());

	auto stone1Tex = mTextures["stone1Tex"]->Resource;
	auto woodTex = mTextures["woodTex"]->Resource;
	auto grassTex = mTextures["grassTex"]->Resource;
	auto roofTex = mTextures["roofTex"]->Resource;
	auto towerTex = mTextures["towerTex"]->Resource;
	auto triprismTex = mTextures["triprismTex"]->Resource;
	auto torusTex = mTextures["torusTex"]->Resource;
	auto diamondTex = mTextures["diamondTex"]->Resource;

	D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc = {};
	srvDesc.Shader4ComponentMapping = D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
	srvDesc.Format = stone1Tex->GetDesc().Format;
	srvDesc.ViewDimension = D3D12_SRV_DIMENSION_TEXTURE2D;
	srvDesc.Texture2D.MostDetailedMip = 0;
	srvDesc.Texture2D.MipLevels = stone1Tex->GetDesc().MipLevels;
	srvDesc.Texture2D.ResourceMinLODClamp = 0.0f;

	// stone
	md3dDevice->CreateShaderResourceView(stone1Tex.Get(), &srvDesc, hDescriptor);

	// wood
	hDescriptor.Offset(1, mCbvSrvUavDescriptorSize);
	srvDesc.Format = woodTex->GetDesc().Format;
	srvDesc.Texture2D.MipLevels = woodTex->GetDesc().MipLevels;
	md3dDevice->CreateShaderResourceView(woodTex.Get(), &srvDesc, hDescriptor);

	// grass
	hDescriptor.Offset(1, mCbvSrvUavDescriptorSize);
	srvDesc.Format = grassTex->GetDesc().Format;
	srvDesc.Texture2D.MipLevels = grassTex->GetDesc().MipLevels;
	md3dDevice->CreateShaderResourceView(grassTex.Get(), &srvDesc, hDescriptor);

	// roof
	hDescriptor.Offset(1, mCbvSrvUavDescriptorSize);
	srvDesc.Format = roofTex->GetDesc().Format;
	srvDesc.Texture2D.MipLevels = roofTex->GetDesc().MipLevels;
	md3dDevice->CreateShaderResourceView(roofTex.Get(), &srvDesc, hDescriptor);

	// tower
	hDescriptor.Offset(1, mCbvSrvUavDescriptorSize);
	srvDesc.Format = towerTex->GetDesc().Format;
	srvDesc.Texture2D.MipLevels = towerTex->GetDesc().MipLevels;
	md3dDevice->CreateShaderResourceView(towerTex.Get(), &srvDesc, hDescriptor);

	// triPrism
	hDescriptor.Offset(1, mCbvSrvUavDescriptorSize);
	srvDesc.Format = triprismTex->GetDesc().Format;
	srvDesc.Texture2D.MipLevels = triprismTex->GetDesc().MipLevels;
	md3dDevice->CreateShaderResourceView(triprismTex.Get(), &srvDesc, hDescriptor);

	// torus
	hDescriptor.Offset(1, mCbvSrvUavDescriptorSize);
	srvDesc.Format = torusTex->GetDesc().Format;
	srvDesc.Texture2D.MipLevels = torusTex->GetDesc().MipLevels;
	md3dDevice->CreateShaderResourceView(torusTex.Get(), &srvDesc, hDescriptor);

	// diamond
	hDescriptor.Offset(1, mCbvSrvUavDescriptorSize);
	srvDesc.Format = diamondTex->GetDesc().Format;
	srvDesc.Texture2D.MipLevels = diamondTex->GetDesc().MipLevels;
	md3dDevice->CreateShaderResourceView(diamondTex.Get(), &srvDesc, hDescriptor);

	auto waterTex = mTextures["waterTex"]->Resource;

	hDescriptor.Offset(1, mCbvSrvUavDescriptorSize);
	srvDesc.Format = waterTex->GetDesc().Format;
	srvDesc.Texture2D.MipLevels = waterTex->GetDesc().MipLevels;
	md3dDevice->CreateShaderResourceView(waterTex.Get(), &srvDesc, hDescriptor);

	// CBV heap for objects + pass
	UINT objCount = (UINT)mOpaqueRitems.size();
	UINT numDescriptors = (objCount + 1) * gNumFrameResources;

	mPassCbvOffset = objCount * gNumFrameResources;

	D3D12_DESCRIPTOR_HEAP_DESC cbvHeapDesc;
	cbvHeapDesc.NumDescriptors = numDescriptors;
	cbvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
	cbvHeapDesc.Flags = D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
	cbvHeapDesc.NodeMask = 0;
	ThrowIfFailed(md3dDevice->CreateDescriptorHeap(&cbvHeapDesc, IID_PPV_ARGS(&mCbvHeap)));
}

void ShapesApp::BuildConstantBufferViews()
{
	UINT objCBByteSize = d3dUtil::CalcConstantBufferByteSize(sizeof(ObjectConstants));

	UINT objCount = (UINT)mOpaqueRitems.size();

	// Need a CBV descriptor for each object for each frame resource.
	for (int frameIndex = 0; frameIndex < gNumFrameResources; ++frameIndex)
	{
		auto objectCB = mFrameResources[frameIndex]->ObjectCB->Resource();
		for (UINT i = 0; i < objCount; ++i)
		{
			D3D12_GPU_VIRTUAL_ADDRESS cbAddress = objectCB->GetGPUVirtualAddress();

			// Offset to the ith object constant buffer in the buffer.
			cbAddress += i * objCBByteSize;

			// Offset to the object cbv in the descriptor heap.
			int heapIndex = frameIndex * objCount + i;
			auto handle = CD3DX12_CPU_DESCRIPTOR_HANDLE(mCbvHeap->GetCPUDescriptorHandleForHeapStart());
			handle.Offset(heapIndex, mCbvSrvUavDescriptorSize);

			D3D12_CONSTANT_BUFFER_VIEW_DESC cbvDesc;
			cbvDesc.BufferLocation = cbAddress;
			cbvDesc.SizeInBytes = objCBByteSize;

			md3dDevice->CreateConstantBufferView(&cbvDesc, handle);
		}
	}

	UINT passCBByteSize = d3dUtil::CalcConstantBufferByteSize(sizeof(PassConstants));

	// Last three descriptors are the pass CBVs for each frame resource.
	for (int frameIndex = 0; frameIndex < gNumFrameResources; ++frameIndex)
	{
		auto passCB = mFrameResources[frameIndex]->PassCB->Resource();
		D3D12_GPU_VIRTUAL_ADDRESS cbAddress = passCB->GetGPUVirtualAddress();

		// Offset to the pass cbv in the descriptor heap.
		int heapIndex = mPassCbvOffset + frameIndex;
		auto handle = CD3DX12_CPU_DESCRIPTOR_HANDLE(mCbvHeap->GetCPUDescriptorHandleForHeapStart());
		handle.Offset(heapIndex, mCbvSrvUavDescriptorSize);

		D3D12_CONSTANT_BUFFER_VIEW_DESC cbvDesc;
		cbvDesc.BufferLocation = cbAddress;
		cbvDesc.SizeInBytes = passCBByteSize;

		md3dDevice->CreateConstantBufferView(&cbvDesc, handle);
	}
}

// NEW: default samplers used by the texture shader
std::array<const CD3DX12_STATIC_SAMPLER_DESC, 6> ShapesApp::GetStaticSamplers()
{
	const CD3DX12_STATIC_SAMPLER_DESC pointWrap(
		0, D3D12_FILTER_MIN_MAG_MIP_POINT,
		D3D12_TEXTURE_ADDRESS_MODE_WRAP,
		D3D12_TEXTURE_ADDRESS_MODE_WRAP,
		D3D12_TEXTURE_ADDRESS_MODE_WRAP);

	const CD3DX12_STATIC_SAMPLER_DESC pointClamp(
		1, D3D12_FILTER_MIN_MAG_MIP_POINT,
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP,
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP,
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP);

	const CD3DX12_STATIC_SAMPLER_DESC linearWrap(
		2, D3D12_FILTER_MIN_MAG_MIP_LINEAR,
		D3D12_TEXTURE_ADDRESS_MODE_WRAP,
		D3D12_TEXTURE_ADDRESS_MODE_WRAP,
		D3D12_TEXTURE_ADDRESS_MODE_WRAP);

	const CD3DX12_STATIC_SAMPLER_DESC linearClamp(
		3, D3D12_FILTER_MIN_MAG_MIP_LINEAR,
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP,
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP,
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP);

	const CD3DX12_STATIC_SAMPLER_DESC anisotropicWrap(
		4, D3D12_FILTER_ANISOTROPIC,
		D3D12_TEXTURE_ADDRESS_MODE_WRAP,
		D3D12_TEXTURE_ADDRESS_MODE_WRAP,
		D3D12_TEXTURE_ADDRESS_MODE_WRAP,
		0.0f, 8);

	const CD3DX12_STATIC_SAMPLER_DESC anisotropicClamp(
		5, D3D12_FILTER_ANISOTROPIC,
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP,
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP,
		D3D12_TEXTURE_ADDRESS_MODE_CLAMP,
		0.0f, 8);

	return { pointWrap, pointClamp, linearWrap, linearClamp, anisotropicWrap, anisotropicClamp };
}

void ShapesApp::BuildRootSignature()
{
	CD3DX12_DESCRIPTOR_RANGE texTable;
	texTable.Init(D3D12_DESCRIPTOR_RANGE_TYPE_SRV, 1, 0);

	CD3DX12_ROOT_PARAMETER slotRootParameter[4];

	slotRootParameter[0].InitAsDescriptorTable(1, &texTable, D3D12_SHADER_VISIBILITY_PIXEL);
	slotRootParameter[1].InitAsConstantBufferView(0);
	slotRootParameter[2].InitAsConstantBufferView(1);
	slotRootParameter[3].InitAsConstantBufferView(2);

	auto staticSamplers = GetStaticSamplers();

	CD3DX12_ROOT_SIGNATURE_DESC rootSigDesc(
		4,
		slotRootParameter,
		(UINT)staticSamplers.size(),
		staticSamplers.data(),
		D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT);

	ComPtr<ID3DBlob> serializedRootSig = nullptr;
	ComPtr<ID3DBlob> errorBlob = nullptr;
	HRESULT hr = D3D12SerializeRootSignature(
		&rootSigDesc,
		D3D_ROOT_SIGNATURE_VERSION_1,
		serializedRootSig.GetAddressOf(),
		errorBlob.GetAddressOf());

	if (errorBlob != nullptr)
		::OutputDebugStringA((char*)errorBlob->GetBufferPointer());

	ThrowIfFailed(hr);

	ThrowIfFailed(md3dDevice->CreateRootSignature(
		0,
		serializedRootSig->GetBufferPointer(),
		serializedRootSig->GetBufferSize(),
		IID_PPV_ARGS(mRootSignature.GetAddressOf())));
}
void ShapesApp::BuildShadersAndInputLayout()
{
	mShaders["standardVS"] = d3dUtil::CompileShader(L"Tex.hlsl", nullptr, "VS", "vs_5_1");
	mShaders["opaquePS"] = d3dUtil::CompileShader(L"Tex.hlsl", nullptr, "PS", "ps_5_1");

	// CHANGED: now using position, normal, and UVs for textures
	mInputLayout =
	{
		{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
		{ "NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 12, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
		{ "TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT, 0, 24, D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
	};
}


// ============================================================
// PART 1: CUSTOM SHAPE GEOMETRY
// This function creates all base and custom shapes used in the
// castle scene, including cone, pyramid, wedge, diamond,
// triangular prism, and torus.
// It also packs all vertices and indices into one large shared
// vertex/index buffer and registers each shape in DrawArgs.
// ============================================================
void ShapesApp::BuildShapeGeometry()
{
	GeometryGenerator geoGen;

	GeometryGenerator::MeshData box = geoGen.CreateBox(1.0f, 1.0f, 1.0f, 0);
	GeometryGenerator::MeshData grid = geoGen.CreateGrid(50.0f, 50.0f, 60, 80);
	GeometryGenerator::MeshData sphere = geoGen.CreateSphere(0.5f, 20, 20);
	GeometryGenerator::MeshData cylinder = geoGen.CreateCylinder(0.5f, 0.3f, 3.0f, 20, 20);

	// NEW primitives (Part 1)
	GeometryGenerator::MeshData cone = geoGen.CreateCone(0.6f, 1.2f, 24, 8);
	GeometryGenerator::MeshData pyramid = geoGen.CreatePyramid(1.0f, 1.2f, 1.0f);
	GeometryGenerator::MeshData wedge = geoGen.CreateWedge(1.2f, 1.0f, 1.2f);
	GeometryGenerator::MeshData diamond = geoGen.CreateDiamond(1.0f, 1.4f, 1.0f);
	GeometryGenerator::MeshData triPrism = geoGen.CreateTriangularPrism(1.2f, 1.0f, 1.2f);
	GeometryGenerator::MeshData torus = geoGen.CreateTorus(0.9f, 0.25f, 32, 16);
	GeometryGenerator::MeshData water = geoGen.CreateGrid(80.0f, 80.0f, 60, 40);

	// ----------------------------
	// Vertex offsets
	UINT boxVertexOffset = 0;
	UINT gridVertexOffset = boxVertexOffset + (UINT)box.Vertices.size();
	UINT sphereVertexOffset = gridVertexOffset + (UINT)grid.Vertices.size();
	UINT cylinderVertexOffset = sphereVertexOffset + (UINT)sphere.Vertices.size();

	UINT coneVertexOffset = cylinderVertexOffset + (UINT)cylinder.Vertices.size();
	UINT pyramidVertexOffset = coneVertexOffset + (UINT)cone.Vertices.size();
	UINT wedgeVertexOffset = pyramidVertexOffset + (UINT)pyramid.Vertices.size();
	UINT diamondVertexOffset = wedgeVertexOffset + (UINT)wedge.Vertices.size();
	UINT triPrismVertexOffset = diamondVertexOffset + (UINT)diamond.Vertices.size();
	UINT torusVertexOffset = triPrismVertexOffset + (UINT)triPrism.Vertices.size();
	UINT waterVertexOffset = torusVertexOffset + (UINT)torus.Vertices.size();

	// Index offsets (use Indices32 sizes  but well pack 16-bit later)
	UINT boxIndexOffset = 0;
	UINT gridIndexOffset = boxIndexOffset + (UINT)box.Indices32.size();
	UINT sphereIndexOffset = gridIndexOffset + (UINT)grid.Indices32.size();
	UINT cylinderIndexOffset = sphereIndexOffset + (UINT)sphere.Indices32.size();

	UINT coneIndexOffset = cylinderIndexOffset + (UINT)cylinder.Indices32.size();
	UINT pyramidIndexOffset = coneIndexOffset + (UINT)cone.Indices32.size();
	UINT wedgeIndexOffset = pyramidIndexOffset + (UINT)pyramid.Indices32.size();
	UINT diamondIndexOffset = wedgeIndexOffset + (UINT)wedge.Indices32.size();
	UINT triPrismIndexOffset = diamondIndexOffset + (UINT)diamond.Indices32.size();
	UINT torusIndexOffset = triPrismIndexOffset + (UINT)triPrism.Indices32.size();
	UINT waterIndexOffset = torusIndexOffset + (UINT)torus.Indices32.size();

	// ----------------------------
	// Submeshes
	SubmeshGeometry boxSubmesh;
	boxSubmesh.IndexCount = (UINT)box.Indices32.size();
	boxSubmesh.StartIndexLocation = boxIndexOffset;
	boxSubmesh.BaseVertexLocation = boxVertexOffset;

	SubmeshGeometry gridSubmesh;
	gridSubmesh.IndexCount = (UINT)grid.Indices32.size();
	gridSubmesh.StartIndexLocation = gridIndexOffset;
	gridSubmesh.BaseVertexLocation = gridVertexOffset;

	SubmeshGeometry sphereSubmesh;
	sphereSubmesh.IndexCount = (UINT)sphere.Indices32.size();
	sphereSubmesh.StartIndexLocation = sphereIndexOffset;
	sphereSubmesh.BaseVertexLocation = sphereVertexOffset;

	SubmeshGeometry cylinderSubmesh;
	cylinderSubmesh.IndexCount = (UINT)cylinder.Indices32.size();
	cylinderSubmesh.StartIndexLocation = cylinderIndexOffset;
	cylinderSubmesh.BaseVertexLocation = cylinderVertexOffset;

	SubmeshGeometry coneSubmesh;
	coneSubmesh.IndexCount = (UINT)cone.Indices32.size();
	coneSubmesh.StartIndexLocation = coneIndexOffset;
	coneSubmesh.BaseVertexLocation = coneVertexOffset;

	SubmeshGeometry pyramidSubmesh;
	pyramidSubmesh.IndexCount = (UINT)pyramid.Indices32.size();
	pyramidSubmesh.StartIndexLocation = pyramidIndexOffset;
	pyramidSubmesh.BaseVertexLocation = pyramidVertexOffset;

	SubmeshGeometry wedgeSubmesh;
	wedgeSubmesh.IndexCount = (UINT)wedge.Indices32.size();
	wedgeSubmesh.StartIndexLocation = wedgeIndexOffset;
	wedgeSubmesh.BaseVertexLocation = wedgeVertexOffset;

	SubmeshGeometry diamondSubmesh;
	diamondSubmesh.IndexCount = (UINT)diamond.Indices32.size();
	diamondSubmesh.StartIndexLocation = diamondIndexOffset;
	diamondSubmesh.BaseVertexLocation = diamondVertexOffset;

	SubmeshGeometry triPrismSubmesh;
	triPrismSubmesh.IndexCount = (UINT)triPrism.Indices32.size();
	triPrismSubmesh.StartIndexLocation = triPrismIndexOffset;
	triPrismSubmesh.BaseVertexLocation = triPrismVertexOffset;

	SubmeshGeometry torusSubmesh;
	torusSubmesh.IndexCount = (UINT)torus.Indices32.size();
	torusSubmesh.StartIndexLocation = torusIndexOffset;
	torusSubmesh.BaseVertexLocation = torusVertexOffset;

	// PART 3: Water plane mesh
	SubmeshGeometry waterSubmesh;
	waterSubmesh.IndexCount = (UINT)water.Indices32.size();
	waterSubmesh.StartIndexLocation = waterIndexOffset;
	waterSubmesh.BaseVertexLocation = waterVertexOffset;


	// ----------------------------
	// Total vertices
	auto totalVertexCount =
		box.Vertices.size() +
		grid.Vertices.size() +
		sphere.Vertices.size() +
		cylinder.Vertices.size() +
		cone.Vertices.size() +
		pyramid.Vertices.size() +
		wedge.Vertices.size() +
		diamond.Vertices.size() +
		triPrism.Vertices.size() +
		torus.Vertices.size() +
		water.Vertices.size();

	std::vector<Vertex> vertices(totalVertexCount);

	UINT k = 0;

	// CHANGED: copy position, normal, and UV instead of color
	auto CopyVerts = [&](const GeometryGenerator::MeshData& src)
		{
			for (size_t i = 0; i < src.Vertices.size(); ++i, ++k)
			{
				vertices[k].Pos = src.Vertices[i].Position;
				vertices[k].Normal = src.Vertices[i].Normal;
				vertices[k].TexC = src.Vertices[i].TexC;
			}
		};

	CopyVerts(box);
	CopyVerts(grid);
	CopyVerts(sphere);
	CopyVerts(cylinder);
	CopyVerts(cone);
	CopyVerts(pyramid);
	CopyVerts(wedge);
	CopyVerts(diamond);
	CopyVerts(triPrism);
	CopyVerts(torus);
	CopyVerts(water);

	// ----------------------------
	// Indices (16-bit)
	std::vector<std::uint16_t> indices;
	indices.insert(indices.end(), std::begin(box.GetIndices16()), std::end(box.GetIndices16()));
	indices.insert(indices.end(), std::begin(grid.GetIndices16()), std::end(grid.GetIndices16()));
	indices.insert(indices.end(), std::begin(sphere.GetIndices16()), std::end(sphere.GetIndices16()));
	indices.insert(indices.end(), std::begin(cylinder.GetIndices16()), std::end(cylinder.GetIndices16()));

	indices.insert(indices.end(), std::begin(cone.GetIndices16()), std::end(cone.GetIndices16()));
	indices.insert(indices.end(), std::begin(pyramid.GetIndices16()), std::end(pyramid.GetIndices16()));
	indices.insert(indices.end(), std::begin(wedge.GetIndices16()), std::end(wedge.GetIndices16()));
	indices.insert(indices.end(), std::begin(diamond.GetIndices16()), std::end(diamond.GetIndices16()));
	indices.insert(indices.end(), std::begin(triPrism.GetIndices16()), std::end(triPrism.GetIndices16()));
	indices.insert(indices.end(), std::begin(torus.GetIndices16()), std::end(torus.GetIndices16()));
	indices.insert(indices.end(), std::begin(water.GetIndices16()), std::end(water.GetIndices16()));

	const UINT vbByteSize = (UINT)vertices.size() * sizeof(Vertex);
	const UINT ibByteSize = (UINT)indices.size() * sizeof(std::uint16_t);

	auto geo = std::make_unique<MeshGeometry>();
	geo->Name = "shapeGeo";

	ThrowIfFailed(D3DCreateBlob(vbByteSize, &geo->VertexBufferCPU));
	CopyMemory(geo->VertexBufferCPU->GetBufferPointer(), vertices.data(), vbByteSize);

	ThrowIfFailed(D3DCreateBlob(ibByteSize, &geo->IndexBufferCPU));
	CopyMemory(geo->IndexBufferCPU->GetBufferPointer(), indices.data(), ibByteSize);

	geo->VertexBufferGPU = d3dUtil::CreateDefaultBuffer(md3dDevice.Get(),
		mCommandList.Get(), vertices.data(), vbByteSize, geo->VertexBufferUploader);

	geo->IndexBufferGPU = d3dUtil::CreateDefaultBuffer(md3dDevice.Get(),
		mCommandList.Get(), indices.data(), ibByteSize, geo->IndexBufferUploader);

	geo->VertexByteStride = sizeof(Vertex);
	geo->VertexBufferByteSize = vbByteSize;
	geo->IndexFormat = DXGI_FORMAT_R16_UINT;
	geo->IndexBufferByteSize = ibByteSize;

	// Register draw args
	geo->DrawArgs["box"] = boxSubmesh;
	geo->DrawArgs["grid"] = gridSubmesh;
	geo->DrawArgs["sphere"] = sphereSubmesh;
	geo->DrawArgs["cylinder"] = cylinderSubmesh;

	// NEW draw args
	geo->DrawArgs["cone"] = coneSubmesh;
	geo->DrawArgs["pyramid"] = pyramidSubmesh;
	geo->DrawArgs["wedge"] = wedgeSubmesh;
	geo->DrawArgs["diamond"] = diamondSubmesh;
	geo->DrawArgs["triPrism"] = triPrismSubmesh;
	geo->DrawArgs["torus"] = torusSubmesh;
	geo->DrawArgs["water"] = waterSubmesh;

	mGeometries[geo->Name] = std::move(geo);
}

// ============================================================
// PART 3: WATER MATERIAL
// Creates scene materials, including the water material.
// The alpha value in DiffuseAlbedo controls how transparent
// the water appears.
// ============================================================
void ShapesApp::BuildMaterials()
{
	auto stone1 = std::make_unique<Material>();
	stone1->Name = "stone1";
	stone1->MatCBIndex = 0;
	stone1->DiffuseSrvHeapIndex = 0;
	stone1->DiffuseAlbedo = XMFLOAT4(1, 1, 1, 1);
	stone1->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	stone1->Roughness = 0.3f;
	mMaterials["stone1"] = std::move(stone1);

	auto wood = std::make_unique<Material>();
	wood->Name = "wood";
	wood->MatCBIndex = 1;
	wood->DiffuseSrvHeapIndex = 1;
	wood->DiffuseAlbedo = XMFLOAT4(1, 1, 1, 1);
	wood->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	wood->Roughness = 0.4f;
	mMaterials["wood"] = std::move(wood);

	auto grass = std::make_unique<Material>();
	grass->Name = "grass";
	grass->MatCBIndex = 2;
	grass->DiffuseSrvHeapIndex = 2;
	grass->DiffuseAlbedo = XMFLOAT4(1, 1, 1, 1);
	grass->FresnelR0 = XMFLOAT3(0.01f, 0.01f, 0.01f);
	grass->Roughness = 0.8f;
	mMaterials["grass"] = std::move(grass);

	auto roof = std::make_unique<Material>();
	roof->Name = "roof";
	roof->MatCBIndex = 3;
	roof->DiffuseSrvHeapIndex = 3;
	roof->DiffuseAlbedo = XMFLOAT4(1, 1, 1, 1);
	roof->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	roof->Roughness = 0.3f;
	mMaterials["roof"] = std::move(roof);

	auto tower = std::make_unique<Material>();
	tower->Name = "tower";
	tower->MatCBIndex = 4;
	tower->DiffuseSrvHeapIndex = 4;
	tower->DiffuseAlbedo = XMFLOAT4(1, 1, 1, 1);
	tower->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	tower->Roughness = 0.6f;
	mMaterials["tower"] = std::move(tower);

	auto triprism = std::make_unique<Material>();
	triprism->Name = "triprism";
	triprism->MatCBIndex = 5;
	triprism->DiffuseSrvHeapIndex = 5;
	triprism->DiffuseAlbedo = XMFLOAT4(1, 1, 1, 1);
	triprism->FresnelR0 = XMFLOAT3(0.02f, 0.02f, 0.02f);
	triprism->Roughness = 0.7f;
	mMaterials["triprism"] = std::move(triprism);

	auto torus = std::make_unique<Material>();
	torus->Name = "torus";
	torus->MatCBIndex = 6;
	torus->DiffuseSrvHeapIndex = 6;
	torus->DiffuseAlbedo = XMFLOAT4(1, 1, 1, 1);
	torus->FresnelR0 = XMFLOAT3(0.08f, 0.06f, 0.02f);
	torus->Roughness = 0.2f;
	mMaterials["torus"] = std::move(torus);

	auto diamond = std::make_unique<Material>();
	diamond->Name = "diamond";
	diamond->MatCBIndex = 7;
	diamond->DiffuseSrvHeapIndex = 7;
	diamond->DiffuseAlbedo = XMFLOAT4(1, 1, 1, 1);
	diamond->FresnelR0 = XMFLOAT3(0.2f, 0.2f, 0.2f);
	diamond->Roughness = 0.05f;
	mMaterials["diamond"] = std::move(diamond);

	auto water = std::make_unique<Material>();
	water->Name = "water";
	water->MatCBIndex = 8;
	water->DiffuseSrvHeapIndex = 8;
	water->DiffuseAlbedo = XMFLOAT4(1.0f, 1.0f, 1.0f, 0.7f);
	water->FresnelR0 = XMFLOAT3(0.1f, 0.1f, 0.1f);
	water->Roughness = 0.0f;
	mMaterials["water"] = std::move(water);
}

// ============================================================
// PART 3: TRANSPARENCY PIPELINE STATE
// Creates the PSOs used for rendering.
// Includes a transparent PSO for the water plane using
// alpha blending (SRC_ALPHA, INV_SRC_ALPHA).
// ============================================================
void ShapesApp::BuildPSOs()
{
	D3D12_GRAPHICS_PIPELINE_STATE_DESC opaquePsoDesc;

	ZeroMemory(&opaquePsoDesc, sizeof(D3D12_GRAPHICS_PIPELINE_STATE_DESC));
	opaquePsoDesc.InputLayout = { mInputLayout.data(), (UINT)mInputLayout.size() };
	opaquePsoDesc.pRootSignature = mRootSignature.Get();

	opaquePsoDesc.VS =
	{
		reinterpret_cast<BYTE*>(mShaders["standardVS"]->GetBufferPointer()),
		mShaders["standardVS"]->GetBufferSize()
	};

	opaquePsoDesc.PS =
	{
		reinterpret_cast<BYTE*>(mShaders["opaquePS"]->GetBufferPointer()),
		mShaders["opaquePS"]->GetBufferSize()
	};

	opaquePsoDesc.RasterizerState = CD3DX12_RASTERIZER_DESC(D3D12_DEFAULT);
	opaquePsoDesc.BlendState = CD3DX12_BLEND_DESC(D3D12_DEFAULT);
	opaquePsoDesc.DepthStencilState = CD3DX12_DEPTH_STENCIL_DESC(D3D12_DEFAULT);
	opaquePsoDesc.SampleMask = UINT_MAX;
	opaquePsoDesc.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
	opaquePsoDesc.NumRenderTargets = 1;
	opaquePsoDesc.RTVFormats[0] = mBackBufferFormat;
	opaquePsoDesc.SampleDesc.Count = m4xMsaaState ? 4 : 1;
	opaquePsoDesc.SampleDesc.Quality = m4xMsaaState ? (m4xMsaaQuality - 1) : 0;
	opaquePsoDesc.DSVFormat = mDepthStencilFormat;

	// OPAQUE
	ThrowIfFailed(md3dDevice->CreateGraphicsPipelineState(
		&opaquePsoDesc, IID_PPV_ARGS(&mPSOs["opaque"])));

	// WIREFRAME
	D3D12_GRAPHICS_PIPELINE_STATE_DESC opaqueWireframePsoDesc = opaquePsoDesc;
	opaqueWireframePsoDesc.RasterizerState.FillMode = D3D12_FILL_MODE_WIREFRAME;

	ThrowIfFailed(md3dDevice->CreateGraphicsPipelineState(
		&opaqueWireframePsoDesc, IID_PPV_ARGS(&mPSOs["opaque_wireframe"])));

	// TRANSPARENT
	D3D12_GRAPHICS_PIPELINE_STATE_DESC transparentPsoDesc = opaquePsoDesc;

	D3D12_RENDER_TARGET_BLEND_DESC transparencyBlendDesc;
	transparencyBlendDesc.BlendEnable = true;
	transparencyBlendDesc.LogicOpEnable = false;

	transparencyBlendDesc.SrcBlend = D3D12_BLEND_SRC_ALPHA;
	transparencyBlendDesc.DestBlend = D3D12_BLEND_INV_SRC_ALPHA;
	transparencyBlendDesc.BlendOp = D3D12_BLEND_OP_ADD;

	transparencyBlendDesc.SrcBlendAlpha = D3D12_BLEND_ONE;
	transparencyBlendDesc.DestBlendAlpha = D3D12_BLEND_ZERO;
	transparencyBlendDesc.BlendOpAlpha = D3D12_BLEND_OP_ADD;

	transparencyBlendDesc.LogicOp = D3D12_LOGIC_OP_NOOP;
	transparencyBlendDesc.RenderTargetWriteMask = D3D12_COLOR_WRITE_ENABLE_ALL;

	transparentPsoDesc.BlendState.RenderTarget[0] = transparencyBlendDesc;

	ThrowIfFailed(md3dDevice->CreateGraphicsPipelineState(
		&transparentPsoDesc, IID_PPV_ARGS(&mPSOs["transparent"])));
}



void ShapesApp::BuildFrameResources()
{
	for (int i = 0; i < gNumFrameResources; ++i)
	{
		// CHANGED: also pass in number of materials
		mFrameResources.push_back(std::make_unique<FrameResource>(md3dDevice.Get(),
			1, (UINT)mAllRitems.size(), (UINT)mMaterials.size()));
	}
}


// ============================================================
// PART 1: RENDER ITEM HELPER
// Helper function used to turn a shape/submesh into a renderable
// object with a world transform, material, and object constant
// buffer index.
// This makes castle construction cleaner inside BuildRenderItems().
// ============================================================
static std::unique_ptr<RenderItem> MakeShapeRitem(
	MeshGeometry* geo,
	Material* mat,
	const std::string& drawArgName,
	UINT objCBIndex,
	const DirectX::XMMATRIX& world)
{
	auto ritem = std::make_unique<RenderItem>();
	XMStoreFloat4x4(&ritem->World, world);

	ritem->ObjCBIndex = objCBIndex;
	ritem->Mat = mat;   // NEW: assign material to this object
	ritem->Geo = geo;
	ritem->PrimitiveType = D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST;
	ritem->IndexCount = geo->DrawArgs[drawArgName].IndexCount;
	ritem->StartIndexLocation = geo->DrawArgs[drawArgName].StartIndexLocation;
	ritem->BaseVertexLocation = geo->DrawArgs[drawArgName].BaseVertexLocation;

	return ritem;
}


// ============================================================
// PART 1: CASTLE CONSTRUCTION
// This function builds the full castle scene by placing the
// shapes created in BuildShapeGeometry() into the world.
// It positions the ground, towers, walls, gate, decorations,
// battlements, and later also adds the water plane.
// ============================================================
void ShapesApp::BuildRenderItems()
{
	MeshGeometry* geo = mGeometries["shapeGeo"].get();
	UINT objCBIndex = 0;

	// CHANGED: now each object also gets a material
	auto Add = [&](const std::string& name, const std::string& matName, const XMMATRIX& world)
		{
			mAllRitems.push_back(
				MakeShapeRitem(geo, mMaterials[matName].get(), name, objCBIndex++, world)
			);
		};

	// Ground
	Add("grid", "grass", XMMatrixTranslation(0.0f, -2.4f, 0.0f));

	// Hill base under castle (gives elevation)
Add("box", "stone1", XMMatrixScaling(28.0f, 4.0f, 28.0f) * XMMatrixTranslation(0.0f, -2.0f, 0.0f));


	// NEW: repeat the grass texture across the ground
	XMStoreFloat4x4(&mAllRitems.back()->TexTransform, XMMatrixScaling(8.0f, 8.0f, 1.0f));

	// === Trees ===
	auto AddTree = [&](float x, float z)
		{
			float trunkHeight = 2.0f;

			// Trunk
			Add("cylinder", "wood",
				XMMatrixScaling(0.3f, trunkHeight, 0.3f) *
				XMMatrixTranslation(x, trunkHeight * 0.5f, z));

			// Leaves
			Add("cone", "grass",
				XMMatrixScaling(1.2f, 2.5f, 1.2f) *
				XMMatrixTranslation(x, trunkHeight + 1.5f, z));
		};

	AddTree(-15.0f, -10.0f);
	AddTree(15.0f, -10.0f);
	AddTree(-20.0f, 5.0f);
	AddTree(20.0f, 5.0f);
	AddTree(0.0f, 20.0f);
	AddTree(-20.0f, -15.0f);
	AddTree(20.0f, -15.0f);
	AddTree(-23.0f, 5.0f);
	AddTree(25.0f, 5.0f);
	AddTree(15.0f, 20.0f);

	// === Towers ===
	float towerHeight = 4.0f;
	float towerRadius = 2.0f;
	float towerY = towerHeight * 0.5f;

	auto Tower = [&](float x, float z)
		{
			Add("cylinder", "tower",
				XMMatrixScaling(towerRadius, towerHeight, towerRadius) *
				XMMatrixTranslation(x, towerY, z));

			Add("cone", "roof",
				XMMatrixScaling(towerRadius * 1.2f, 4.0f, towerRadius * 1.2f) *
				XMMatrixTranslation(x, towerHeight + 1.2f, z));

			Add("diamond", "diamond",
				XMMatrixScaling(1.2f, 1.2f, 1.2f) *
				XMMatrixTranslation(x, towerHeight + 7.0f, z));
		};

	Tower(-8, -8);
	Tower(8, -8);
	Tower(8, 8);
	Tower(-8, 8);

	// === Walls ===
	float wallHeight = 3.0f;
	float wallY = wallHeight * 0.5f;

	float wallLen = 16.0f;
	float gateWidth = 5.0f;

	float sideLen = (wallLen - gateWidth) * 0.5f;
	float sideOffset = (gateWidth * 0.5f) + (sideLen * 0.5f);

	float battH = 1.0f;
	float battY = wallHeight + battH * 0.5f;

	float start = -6.0f;
	float end = 6.0f;
	float step = 3.0f;

	Add("box", "stone1",
		XMMatrixScaling(sideLen, wallHeight, 1) *
		XMMatrixTranslation(-sideOffset, wallY, -8));

	Add("box", "stone1",
		XMMatrixScaling(sideLen, wallHeight, 1) *
		XMMatrixTranslation(sideOffset, wallY, -8));

	Add("box", "stone1",
		XMMatrixScaling(16, wallHeight, 1) *
		XMMatrixTranslation(0, wallY, 8));

	Add("box", "stone1",
		XMMatrixScaling(16, wallHeight, 1) *
		XMMatrixRotationY(XM_PIDIV2) *
		XMMatrixTranslation(-8, wallY, 0));

	Add("box", "stone1",
		XMMatrixScaling(16, wallHeight, 1) *
		XMMatrixRotationY(XM_PIDIV2) *
		XMMatrixTranslation(8, wallY, 0));

	// === Gate ===
	Add("wedge", "wood",
		XMMatrixScaling(4, 3, 2) *
	    XMMatrixTranslation(0, 0.5f, -8.2f));

	Add("wedge", "wood",
		XMMatrixScaling(4, 3, 2) *
		XMMatrixRotationY(XM_PI) *
		XMMatrixTranslation(0, 0.5f, -8.2f));

	// Bridge leading to gate
	for (int i = 0; i < 12; i++)
	{
		Add("box", "stone1",
			XMMatrixScaling(3.0f, 0.5f, 2.0f) *
			XMMatrixTranslation(0.0f, 0.5f, -9.0f - i * 2.0f));
	}

	// === Decorations ===
	Add("diamond", "diamond",
		XMMatrixScaling(1.5f, 1.5f, 1.5f) *
		XMMatrixTranslation(0, 7.5f, 0));

	Add("torus", "torus",
		XMMatrixScaling(1.3f, 1.3f, 1.3f) *
		XMMatrixTranslation(0, 5.0f, -9.5f));

	// Battlements
	float gateClearHalfWidth = gateWidth * 0.5f;

	for (float x = start; x <= end; x += step)
	{
		if (fabsf(x) < gateClearHalfWidth) continue;

		Add("triPrism", "triprism",
			XMMatrixScaling(1, battH, 1) *
			XMMatrixTranslation(x, wallHeight, -8.0f));
	}

	for (float x = start; x <= end; x += step)
	{
		Add("triPrism", "triprism",
			XMMatrixScaling(1, battH, 1) *
			XMMatrixTranslation(x, wallHeight, 8.0f));
	}

	for (float z = start; z <= end; z += step)
	{
		Add("triPrism", "triprism",
			XMMatrixScaling(1, battH, 1) *
			XMMatrixRotationY(XM_PIDIV2) *
			XMMatrixTranslation(-8.0f, wallHeight, z));
	}

	for (float z = start; z <= end; z += step)
	{
		Add("triPrism", "triprism",
			XMMatrixScaling(1, battH, 1) *
			XMMatrixRotationY(XM_PIDIV2) *
			XMMatrixTranslation(8.0f, wallHeight, z));
	}

	// ============================================================
	// PART 3: WATER RENDER ITEM
	// Creates the water plane using the grid mesh and assigns the
	// water material. The water is placed slightly below the ground
	// and added to the transparent render list so it is drawn after
	// opaque objects for correct blending.
	// ============================================================
	mAllRitems.push_back(
		MakeShapeRitem(
			geo,
			mMaterials["water"].get(),
			"water",
			objCBIndex++,
			XMMatrixTranslation(0.0f, -2.5f, 0.0f)   // water height
		)
	);

	XMStoreFloat4x4(&mAllRitems.back()->TexTransform, XMMatrixScaling(6.0f, 6.0f, 1.0f));
	mTransparentRitems.push_back(mAllRitems.back().get());


	// Finalize
	for (auto& e : mAllRitems)
	{
		if (e->Mat == mMaterials["water"].get())
			continue;

		mOpaqueRitems.push_back(e.get());
	}
}


void ShapesApp::DrawRenderItems(ID3D12GraphicsCommandList* cmdList, const std::vector<RenderItem*>& ritems)
{
	UINT objCBByteSize = d3dUtil::CalcConstantBufferByteSize(sizeof(ObjectConstants));
	UINT matCBByteSize = d3dUtil::CalcConstantBufferByteSize(sizeof(MaterialConstants));

	auto objectCB = mCurrFrameResource->ObjectCB->Resource();
	auto matCB = mCurrFrameResource->MaterialCB->Resource();

	for (size_t i = 0; i < ritems.size(); ++i)
	{
		auto ri = ritems[i];

		cmdList->IASetVertexBuffers(0, 1, &ri->Geo->VertexBufferView());
		cmdList->IASetIndexBuffer(&ri->Geo->IndexBufferView());
		cmdList->IASetPrimitiveTopology(ri->PrimitiveType);

		// object CBV
		D3D12_GPU_VIRTUAL_ADDRESS objCBAddress = objectCB->GetGPUVirtualAddress();
		objCBAddress += ri->ObjCBIndex * objCBByteSize;

		// texture SRV
		auto tex = CD3DX12_GPU_DESCRIPTOR_HANDLE(mSrvDescriptorHeap->GetGPUDescriptorHandleForHeapStart());
		tex.Offset(ri->Mat->DiffuseSrvHeapIndex, mCbvSrvUavDescriptorSize);

		// material CB address
		D3D12_GPU_VIRTUAL_ADDRESS matCBAddress = matCB->GetGPUVirtualAddress();
		matCBAddress += ri->Mat->MatCBIndex * matCBByteSize;

		// slot 0 = texture
		cmdList->SetGraphicsRootDescriptorTable(0, tex);

		// slot 1 = object CB
		cmdList->SetGraphicsRootConstantBufferView(1, objCBAddress);

		// slot 3 = material CB
		cmdList->SetGraphicsRootConstantBufferView(3, matCBAddress);

		cmdList->DrawIndexedInstanced(ri->IndexCount, 1, ri->StartIndexLocation, ri->BaseVertexLocation, 0);
	}
}