//***************************************************************************************
// GeometryGenerator.cpp by Frank Luna (C) 2011 All Rights Reserved.
//***************************************************************************************

#include "GeometryGenerator.h"
#include <algorithm>
#include <DirectXMath.h>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace DirectX;

GeometryGenerator::MeshData GeometryGenerator::CreateBox(float width, float height, float depth, uint32 numSubdivisions)
{
    MeshData meshData;

    //
	// Create the vertices.
	//

	Vertex v[24];

	float w2 = 0.5f*width;
	float h2 = 0.5f*height;
	float d2 = 0.5f*depth;
    
	// Fill in the front face vertex data.
	v[0] = Vertex(-w2, -h2, -d2, 0.0f, 0.0f, -1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f);
	v[1] = Vertex(-w2, +h2, -d2, 0.0f, 0.0f, -1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f);
	v[2] = Vertex(+w2, +h2, -d2, 0.0f, 0.0f, -1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f);
	v[3] = Vertex(+w2, -h2, -d2, 0.0f, 0.0f, -1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f);

	// Fill in the back face vertex data.
	v[4] = Vertex(-w2, -h2, +d2, 0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f, 1.0f, 1.0f);
	v[5] = Vertex(+w2, -h2, +d2, 0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f, 0.0f, 1.0f);
	v[6] = Vertex(+w2, +h2, +d2, 0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f);
	v[7] = Vertex(-w2, +h2, +d2, 0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f, 1.0f, 0.0f);

	// Fill in the top face vertex data.
	v[8]  = Vertex(-w2, +h2, -d2, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f);
	v[9]  = Vertex(-w2, +h2, +d2, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f);
	v[10] = Vertex(+w2, +h2, +d2, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f);
	v[11] = Vertex(+w2, +h2, -d2, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f);

	// Fill in the bottom face vertex data.
	v[12] = Vertex(-w2, -h2, -d2, 0.0f, -1.0f, 0.0f, -1.0f, 0.0f, 0.0f, 1.0f, 1.0f);
	v[13] = Vertex(+w2, -h2, -d2, 0.0f, -1.0f, 0.0f, -1.0f, 0.0f, 0.0f, 0.0f, 1.0f);
	v[14] = Vertex(+w2, -h2, +d2, 0.0f, -1.0f, 0.0f, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f);
	v[15] = Vertex(-w2, -h2, +d2, 0.0f, -1.0f, 0.0f, -1.0f, 0.0f, 0.0f, 1.0f, 0.0f);

	// Fill in the left face vertex data.
	v[16] = Vertex(-w2, -h2, +d2, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 1.0f);
	v[17] = Vertex(-w2, +h2, +d2, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f);
	v[18] = Vertex(-w2, +h2, -d2, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 1.0f, 0.0f);
	v[19] = Vertex(-w2, -h2, -d2, -1.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 1.0f, 1.0f);

	// Fill in the right face vertex data.
	v[20] = Vertex(+w2, -h2, -d2, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f);
	v[21] = Vertex(+w2, +h2, -d2, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f);
	v[22] = Vertex(+w2, +h2, +d2, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f);
	v[23] = Vertex(+w2, -h2, +d2, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f);

	meshData.Vertices.assign(&v[0], &v[24]);
 
	//
	// Create the indices.
	//

	uint32 i[36];

	// Fill in the front face index data
	i[0] = 0; i[1] = 1; i[2] = 2;
	i[3] = 0; i[4] = 2; i[5] = 3;

	// Fill in the back face index data
	i[6] = 4; i[7]  = 5; i[8]  = 6;
	i[9] = 4; i[10] = 6; i[11] = 7;

	// Fill in the top face index data
	i[12] = 8; i[13] =  9; i[14] = 10;
	i[15] = 8; i[16] = 10; i[17] = 11;

	// Fill in the bottom face index data
	i[18] = 12; i[19] = 13; i[20] = 14;
	i[21] = 12; i[22] = 14; i[23] = 15;

	// Fill in the left face index data
	i[24] = 16; i[25] = 17; i[26] = 18;
	i[27] = 16; i[28] = 18; i[29] = 19;

	// Fill in the right face index data
	i[30] = 20; i[31] = 21; i[32] = 22;
	i[33] = 20; i[34] = 22; i[35] = 23;

	meshData.Indices32.assign(&i[0], &i[36]);

    // Put a cap on the number of subdivisions.
    numSubdivisions = std::min<uint32>(numSubdivisions, 6u);

    for(uint32 i = 0; i < numSubdivisions; ++i)
        Subdivide(meshData);

    return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateSphere(float radius, uint32 sliceCount, uint32 stackCount)
{
    MeshData meshData;

	//
	// Compute the vertices stating at the top pole and moving down the stacks.
	//

	// Poles: note that there will be texture coordinate distortion as there is
	// not a unique point on the texture map to assign to the pole when mapping
	// a rectangular texture onto a sphere.
	Vertex topVertex(0.0f, +radius, 0.0f, 0.0f, +1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f);
	Vertex bottomVertex(0.0f, -radius, 0.0f, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f);

	meshData.Vertices.push_back( topVertex );

	float phiStep   = XM_PI/stackCount;
	float thetaStep = 2.0f*XM_PI/sliceCount;

	// Compute vertices for each stack ring (do not count the poles as rings).
	for(uint32 i = 1; i <= stackCount-1; ++i)
	{
		float phi = i*phiStep;

		// Vertices of ring.
        for(uint32 j = 0; j <= sliceCount; ++j)
		{
			float theta = j*thetaStep;

			Vertex v;

			// spherical to cartesian
			v.Position.x = radius*sinf(phi)*cosf(theta);
			v.Position.y = radius*cosf(phi);
			v.Position.z = radius*sinf(phi)*sinf(theta);

			// Partial derivative of P with respect to theta
			v.TangentU.x = -radius*sinf(phi)*sinf(theta);
			v.TangentU.y = 0.0f;
			v.TangentU.z = +radius*sinf(phi)*cosf(theta);

			XMVECTOR T = XMLoadFloat3(&v.TangentU);
			XMStoreFloat3(&v.TangentU, XMVector3Normalize(T));

			XMVECTOR p = XMLoadFloat3(&v.Position);
			XMStoreFloat3(&v.Normal, XMVector3Normalize(p));

			v.TexC.x = theta / XM_2PI;
			v.TexC.y = phi / XM_PI;

			meshData.Vertices.push_back( v );
		}
	}

	meshData.Vertices.push_back( bottomVertex );

	//
	// Compute indices for top stack.  The top stack was written first to the vertex buffer
	// and connects the top pole to the first ring.
	//

    for(uint32 i = 1; i <= sliceCount; ++i)
	{
		meshData.Indices32.push_back(0);
		meshData.Indices32.push_back(i+1);
		meshData.Indices32.push_back(i);
	}
	
	//
	// Compute indices for inner stacks (not connected to poles).
	//

	// Offset the indices to the index of the first vertex in the first ring.
	// This is just skipping the top pole vertex.
    uint32 baseIndex = 1;
    uint32 ringVertexCount = sliceCount + 1;
	for(uint32 i = 0; i < stackCount-2; ++i)
	{
		for(uint32 j = 0; j < sliceCount; ++j)
		{
			meshData.Indices32.push_back(baseIndex + i*ringVertexCount + j);
			meshData.Indices32.push_back(baseIndex + i*ringVertexCount + j+1);
			meshData.Indices32.push_back(baseIndex + (i+1)*ringVertexCount + j);

			meshData.Indices32.push_back(baseIndex + (i+1)*ringVertexCount + j);
			meshData.Indices32.push_back(baseIndex + i*ringVertexCount + j+1);
			meshData.Indices32.push_back(baseIndex + (i+1)*ringVertexCount + j+1);
		}
	}

	//
	// Compute indices for bottom stack.  The bottom stack was written last to the vertex buffer
	// and connects the bottom pole to the bottom ring.
	//

	// South pole vertex was added last.
	uint32 southPoleIndex = (uint32)meshData.Vertices.size()-1;

	// Offset the indices to the index of the first vertex in the last ring.
	baseIndex = southPoleIndex - ringVertexCount;
	
	for(uint32 i = 0; i < sliceCount; ++i)
	{
		meshData.Indices32.push_back(southPoleIndex);
		meshData.Indices32.push_back(baseIndex+i);
		meshData.Indices32.push_back(baseIndex+i+1);
	}

    return meshData;
}
 
void GeometryGenerator::Subdivide(MeshData& meshData)
{
	// Save a copy of the input geometry.
	MeshData inputCopy = meshData;


	meshData.Vertices.resize(0);
	meshData.Indices32.resize(0);

	//       v1
	//       *
	//      / \
	//     /   \
	//  m0*-----*m1
	//   / \   / \
	//  /   \ /   \
	// *-----*-----*
	// v0    m2     v2

	uint32 numTris = (uint32)inputCopy.Indices32.size()/3;
	for(uint32 i = 0; i < numTris; ++i)
	{
		Vertex v0 = inputCopy.Vertices[ inputCopy.Indices32[i*3+0] ];
		Vertex v1 = inputCopy.Vertices[ inputCopy.Indices32[i*3+1] ];
		Vertex v2 = inputCopy.Vertices[ inputCopy.Indices32[i*3+2] ];

		//
		// Generate the midpoints.
		//

        Vertex m0 = MidPoint(v0, v1);
        Vertex m1 = MidPoint(v1, v2);
        Vertex m2 = MidPoint(v0, v2);

		//
		// Add new geometry.
		//

		meshData.Vertices.push_back(v0); // 0
		meshData.Vertices.push_back(v1); // 1
		meshData.Vertices.push_back(v2); // 2
		meshData.Vertices.push_back(m0); // 3
		meshData.Vertices.push_back(m1); // 4
		meshData.Vertices.push_back(m2); // 5
 
		meshData.Indices32.push_back(i*6+0);
		meshData.Indices32.push_back(i*6+3);
		meshData.Indices32.push_back(i*6+5);

		meshData.Indices32.push_back(i*6+3);
		meshData.Indices32.push_back(i*6+4);
		meshData.Indices32.push_back(i*6+5);

		meshData.Indices32.push_back(i*6+5);
		meshData.Indices32.push_back(i*6+4);
		meshData.Indices32.push_back(i*6+2);

		meshData.Indices32.push_back(i*6+3);
		meshData.Indices32.push_back(i*6+1);
		meshData.Indices32.push_back(i*6+4);
	}
}

GeometryGenerator::Vertex GeometryGenerator::MidPoint(const Vertex& v0, const Vertex& v1)
{
    XMVECTOR p0 = XMLoadFloat3(&v0.Position);
    XMVECTOR p1 = XMLoadFloat3(&v1.Position);

    XMVECTOR n0 = XMLoadFloat3(&v0.Normal);
    XMVECTOR n1 = XMLoadFloat3(&v1.Normal);

    XMVECTOR tan0 = XMLoadFloat3(&v0.TangentU);
    XMVECTOR tan1 = XMLoadFloat3(&v1.TangentU);

    XMVECTOR tex0 = XMLoadFloat2(&v0.TexC);
    XMVECTOR tex1 = XMLoadFloat2(&v1.TexC);

    // Compute the midpoints of all the attributes.  Vectors need to be normalized
    // since linear interpolating can make them not unit length.  
    XMVECTOR pos = 0.5f*(p0 + p1);
    XMVECTOR normal = XMVector3Normalize(0.5f*(n0 + n1));
    XMVECTOR tangent = XMVector3Normalize(0.5f*(tan0+tan1));
    XMVECTOR tex = 0.5f*(tex0 + tex1);

    Vertex v;
    XMStoreFloat3(&v.Position, pos);
    XMStoreFloat3(&v.Normal, normal);
    XMStoreFloat3(&v.TangentU, tangent);
    XMStoreFloat2(&v.TexC, tex);

    return v;
}

GeometryGenerator::MeshData GeometryGenerator::CreateGeosphere(float radius, uint32 numSubdivisions)
{
    MeshData meshData;

	// Put a cap on the number of subdivisions.
    numSubdivisions = std::min<uint32>(numSubdivisions, 6u);

	// Approximate a sphere by tessellating an icosahedron.

	const float X = 0.525731f; 
	const float Z = 0.850651f;

	XMFLOAT3 pos[12] = 
	{
		XMFLOAT3(-X, 0.0f, Z),  XMFLOAT3(X, 0.0f, Z),  
		XMFLOAT3(-X, 0.0f, -Z), XMFLOAT3(X, 0.0f, -Z),    
		XMFLOAT3(0.0f, Z, X),   XMFLOAT3(0.0f, Z, -X), 
		XMFLOAT3(0.0f, -Z, X),  XMFLOAT3(0.0f, -Z, -X),    
		XMFLOAT3(Z, X, 0.0f),   XMFLOAT3(-Z, X, 0.0f), 
		XMFLOAT3(Z, -X, 0.0f),  XMFLOAT3(-Z, -X, 0.0f)
	};

    uint32 k[60] =
	{
		1,4,0,  4,9,0,  4,5,9,  8,5,4,  1,8,4,    
		1,10,8, 10,3,8, 8,3,5,  3,2,5,  3,7,2,    
		3,10,7, 10,6,7, 6,11,7, 6,0,11, 6,1,0, 
		10,1,6, 11,0,9, 2,11,9, 5,2,9,  11,2,7 
	};

    meshData.Vertices.resize(12);
    meshData.Indices32.assign(&k[0], &k[60]);

	for(uint32 i = 0; i < 12; ++i)
		meshData.Vertices[i].Position = pos[i];

	for(uint32 i = 0; i < numSubdivisions; ++i)
		Subdivide(meshData);

	// Project vertices onto sphere and scale.
	for(uint32 i = 0; i < meshData.Vertices.size(); ++i)
	{
		// Project onto unit sphere.
		XMVECTOR n = XMVector3Normalize(XMLoadFloat3(&meshData.Vertices[i].Position));

		// Project onto sphere.
		XMVECTOR p = radius*n;

		XMStoreFloat3(&meshData.Vertices[i].Position, p);
		XMStoreFloat3(&meshData.Vertices[i].Normal, n);

		// Derive texture coordinates from spherical coordinates.
        float theta = atan2f(meshData.Vertices[i].Position.z, meshData.Vertices[i].Position.x);

        // Put in [0, 2pi].
        if(theta < 0.0f)
            theta += XM_2PI;

		float phi = acosf(meshData.Vertices[i].Position.y / radius);

		meshData.Vertices[i].TexC.x = theta/XM_2PI;
		meshData.Vertices[i].TexC.y = phi/XM_PI;

		// Partial derivative of P with respect to theta
		meshData.Vertices[i].TangentU.x = -radius*sinf(phi)*sinf(theta);
		meshData.Vertices[i].TangentU.y = 0.0f;
		meshData.Vertices[i].TangentU.z = +radius*sinf(phi)*cosf(theta);

		XMVECTOR T = XMLoadFloat3(&meshData.Vertices[i].TangentU);
		XMStoreFloat3(&meshData.Vertices[i].TangentU, XMVector3Normalize(T));
	}

    return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateCylinder(float bottomRadius, float topRadius, float height, uint32 sliceCount, uint32 stackCount)
{
    MeshData meshData;

	//
	// Build Stacks.
	// 

	float stackHeight = height / stackCount;

	// Amount to increment radius as we move up each stack level from bottom to top.
	float radiusStep = (topRadius - bottomRadius) / stackCount;

	uint32 ringCount = stackCount+1;

	// Compute vertices for each stack ring starting at the bottom and moving up.
	for(uint32 i = 0; i < ringCount; ++i)
	{
		float y = -0.5f*height + i*stackHeight;
		float r = bottomRadius + i*radiusStep;

		// vertices of ring
		float dTheta = 2.0f*XM_PI/sliceCount;
		for(uint32 j = 0; j <= sliceCount; ++j)
		{
			Vertex vertex;

			float c = cosf(j*dTheta);
			float s = sinf(j*dTheta);

			vertex.Position = XMFLOAT3(r*c, y, r*s);

			vertex.TexC.x = (float)j/sliceCount;
			vertex.TexC.y = 1.0f - (float)i/stackCount;

			// Cylinder can be parameterized as follows, where we introduce v
			// parameter that goes in the same direction as the v tex-coord
			// so that the bitangent goes in the same direction as the v tex-coord.
			//   Let r0 be the bottom radius and let r1 be the top radius.
			//   y(v) = h - hv for v in [0,1].
			//   r(v) = r1 + (r0-r1)v
			//
			//   x(t, v) = r(v)*cos(t)
			//   y(t, v) = h - hv
			//   z(t, v) = r(v)*sin(t)
			// 
			//  dx/dt = -r(v)*sin(t)
			//  dy/dt = 0
			//  dz/dt = +r(v)*cos(t)
			//
			//  dx/dv = (r0-r1)*cos(t)
			//  dy/dv = -h
			//  dz/dv = (r0-r1)*sin(t)

			// This is unit length.
			vertex.TangentU = XMFLOAT3(-s, 0.0f, c);

			float dr = bottomRadius-topRadius;
			XMFLOAT3 bitangent(dr*c, -height, dr*s);

			XMVECTOR T = XMLoadFloat3(&vertex.TangentU);
			XMVECTOR B = XMLoadFloat3(&bitangent);
			XMVECTOR N = XMVector3Normalize(XMVector3Cross(T, B));
			XMStoreFloat3(&vertex.Normal, N);

			meshData.Vertices.push_back(vertex);
		}
	}

	// Add one because we duplicate the first and last vertex per ring
	// since the texture coordinates are different.
	uint32 ringVertexCount = sliceCount+1;

	// Compute indices for each stack.
	for(uint32 i = 0; i < stackCount; ++i)
	{
		for(uint32 j = 0; j < sliceCount; ++j)
		{
			meshData.Indices32.push_back(i*ringVertexCount + j);
			meshData.Indices32.push_back((i+1)*ringVertexCount + j);
			meshData.Indices32.push_back((i+1)*ringVertexCount + j+1);

			meshData.Indices32.push_back(i*ringVertexCount + j);
			meshData.Indices32.push_back((i+1)*ringVertexCount + j+1);
			meshData.Indices32.push_back(i*ringVertexCount + j+1);
		}
	}

	BuildCylinderTopCap(bottomRadius, topRadius, height, sliceCount, stackCount, meshData);
	BuildCylinderBottomCap(bottomRadius, topRadius, height, sliceCount, stackCount, meshData);

    return meshData;
}

void GeometryGenerator::BuildCylinderTopCap(float bottomRadius, float topRadius, float height,
											uint32 sliceCount, uint32 stackCount, MeshData& meshData)
{
	uint32 baseIndex = (uint32)meshData.Vertices.size();

	float y = 0.5f*height;
	float dTheta = 2.0f*XM_PI/sliceCount;

	// Duplicate cap ring vertices because the texture coordinates and normals differ.
	for(uint32 i = 0; i <= sliceCount; ++i)
	{
		float x = topRadius*cosf(i*dTheta);
		float z = topRadius*sinf(i*dTheta);

		// Scale down by the height to try and make top cap texture coord area
		// proportional to base.
		float u = x/height + 0.5f;
		float v = z/height + 0.5f;

		meshData.Vertices.push_back( Vertex(x, y, z, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, u, v) );
	}

	// Cap center vertex.
	meshData.Vertices.push_back( Vertex(0.0f, y, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.5f, 0.5f) );

	// Index of center vertex.
	uint32 centerIndex = (uint32)meshData.Vertices.size()-1;

	for(uint32 i = 0; i < sliceCount; ++i)
	{
		meshData.Indices32.push_back(centerIndex);
		meshData.Indices32.push_back(baseIndex + i+1);
		meshData.Indices32.push_back(baseIndex + i);
	}
}

void GeometryGenerator::BuildCylinderBottomCap(float bottomRadius, float topRadius, float height,
											   uint32 sliceCount, uint32 stackCount, MeshData& meshData)
{
	// 
	// Build bottom cap.
	//

	uint32 baseIndex = (uint32)meshData.Vertices.size();
	float y = -0.5f*height;

	// vertices of ring
	float dTheta = 2.0f*XM_PI/sliceCount;
	for(uint32 i = 0; i <= sliceCount; ++i)
	{
		float x = bottomRadius*cosf(i*dTheta);
		float z = bottomRadius*sinf(i*dTheta);

		// Scale down by the height to try and make top cap texture coord area
		// proportional to base.
		float u = x/height + 0.5f;
		float v = z/height + 0.5f;

		meshData.Vertices.push_back( Vertex(x, y, z, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f, 0.0f, u, v) );
	}

	// Cap center vertex.
	meshData.Vertices.push_back( Vertex(0.0f, y, 0.0f, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.5f, 0.5f) );

	// Cache the index of center vertex.
	uint32 centerIndex = (uint32)meshData.Vertices.size()-1;

	for(uint32 i = 0; i < sliceCount; ++i)
	{
		meshData.Indices32.push_back(centerIndex);
		meshData.Indices32.push_back(baseIndex + i);
		meshData.Indices32.push_back(baseIndex + i+1);
	}
}

GeometryGenerator::MeshData GeometryGenerator::CreateGrid(float width, float depth, uint32 m, uint32 n)
{
    MeshData meshData;

	uint32 vertexCount = m*n;
	uint32 faceCount   = (m-1)*(n-1)*2;

	//
	// Create the vertices.
	//

	float halfWidth = 0.5f*width;
	float halfDepth = 0.5f*depth;

	float dx = width / (n-1);
	float dz = depth / (m-1);

	float du = 1.0f / (n-1);
	float dv = 1.0f / (m-1);

	meshData.Vertices.resize(vertexCount);
	for(uint32 i = 0; i < m; ++i)
	{
		float z = halfDepth - i*dz;
		for(uint32 j = 0; j < n; ++j)
		{
			float x = -halfWidth + j*dx;

			meshData.Vertices[i*n+j].Position = XMFLOAT3(x, 0.0f, z);
			meshData.Vertices[i*n+j].Normal   = XMFLOAT3(0.0f, 1.0f, 0.0f);
			meshData.Vertices[i*n+j].TangentU = XMFLOAT3(1.0f, 0.0f, 0.0f);

			// Stretch texture over grid.
			meshData.Vertices[i*n+j].TexC.x = j*du;
			meshData.Vertices[i*n+j].TexC.y = i*dv;
		}
	}
 
    //
	// Create the indices.
	//

	meshData.Indices32.resize(faceCount*3); // 3 indices per face

	// Iterate over each quad and compute indices.
	uint32 k = 0;
	for(uint32 i = 0; i < m-1; ++i)
	{
		for(uint32 j = 0; j < n-1; ++j)
		{
			meshData.Indices32[k]   = i*n+j;
			meshData.Indices32[k+1] = i*n+j+1;
			meshData.Indices32[k+2] = (i+1)*n+j;

			meshData.Indices32[k+3] = (i+1)*n+j;
			meshData.Indices32[k+4] = i*n+j+1;
			meshData.Indices32[k+5] = (i+1)*n+j+1;

			k += 6; // next quad
		}
	}

    return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateQuad(float x, float y, float w, float h, float depth)
{
    MeshData meshData;

	meshData.Vertices.resize(4);
	meshData.Indices32.resize(6);

	// Position coordinates specified in NDC space.
	meshData.Vertices[0] = Vertex(
        x, y - h, depth,
		0.0f, 0.0f, -1.0f,
		1.0f, 0.0f, 0.0f,
		0.0f, 1.0f);

	meshData.Vertices[1] = Vertex(
		x, y, depth,
		0.0f, 0.0f, -1.0f,
		1.0f, 0.0f, 0.0f,
		0.0f, 0.0f);

	meshData.Vertices[2] = Vertex(
		x+w, y, depth,
		0.0f, 0.0f, -1.0f,
		1.0f, 0.0f, 0.0f,
		1.0f, 0.0f);

	meshData.Vertices[3] = Vertex(
		x+w, y-h, depth,
		0.0f, 0.0f, -1.0f,
		1.0f, 0.0f, 0.0f,
		1.0f, 1.0f);

	meshData.Indices32[0] = 0;
	meshData.Indices32[1] = 1;
	meshData.Indices32[2] = 2;

	meshData.Indices32[3] = 0;
	meshData.Indices32[4] = 2;
	meshData.Indices32[5] = 3;

    return meshData;
}

// Tyron was here |
//				  |
//			      |
//               \ /
//                v
// =======================================================================================
// NEW PRIMITIVES
// this section is the full mesh generation code for new primitives
// =======================================================================================

GeometryGenerator::MeshData GeometryGenerator::CreateCone(float radius, float height, uint32 sliceCount, uint32 stackCount)
{
	MeshData meshData;

	// Side vertices: stacks from bottom (y=0) to top (y=height).
	for (uint32 i = 0; i <= stackCount; ++i)
	{
		float t = (float)i / (float)stackCount;
		float y = t * height;
		float r = radius * (1.0f - t);

		float dTheta = XM_2PI / sliceCount;

		for (uint32 j = 0; j <= sliceCount; ++j)
		{
			float theta = j * dTheta;

			float x = r * cosf(theta);
			float z = r * sinf(theta);

			// Approx cone normal
			XMVECTOR n = XMVectorSet(x, radius / height, z, 0.0f);
			n = XMVector3Normalize(n);

			XMFLOAT3 nf;
			XMStoreFloat3(&nf, n);

			// Tangent around theta
			XMFLOAT3 tangent(-sinf(theta), 0.0f, cosf(theta));

			float u = (float)j / (float)sliceCount;
			float v = 1.0f - t;

			meshData.Vertices.push_back(
				Vertex(x, y, z,
					nf.x, nf.y, nf.z,
					tangent.x, tangent.y, tangent.z,
					u, v)
			);
		}
	}

	// Indices for side
	uint32 ringVerts = sliceCount + 1;
	for (uint32 i = 0; i < stackCount; ++i)
	{
		for (uint32 j = 0; j < sliceCount; ++j)
		{
			uint32 a = i * ringVerts + j;
			uint32 b = (i + 1) * ringVerts + j;
			uint32 c = (i + 1) * ringVerts + (j + 1);
			uint32 d = i * ringVerts + (j + 1);

			meshData.Indices32.push_back(a);
			meshData.Indices32.push_back(b);
			meshData.Indices32.push_back(c);

			meshData.Indices32.push_back(a);
			meshData.Indices32.push_back(c);
			meshData.Indices32.push_back(d);
		}
	}

	// Base cap (y=0)
	uint32 baseCenterIndex = (uint32)meshData.Vertices.size();
	meshData.Vertices.push_back(Vertex(
		0.0f, 0.0f, 0.0f,
		0.0f, -1.0f, 0.0f,
		1.0f, 0.0f, 0.0f,
		0.5f, 0.5f));

	uint32 baseStartIndex = (uint32)meshData.Vertices.size();
	float dTheta = XM_2PI / sliceCount;

	for (uint32 i = 0; i <= sliceCount; ++i)
	{
		float theta = i * dTheta;
		float x = radius * cosf(theta);
		float z = radius * sinf(theta);

		float u = 0.5f + x / (2.0f * radius);
		float v = 0.5f - z / (2.0f * radius);

		meshData.Vertices.push_back(Vertex(
			x, 0.0f, z,
			0.0f, -1.0f, 0.0f,
			1.0f, 0.0f, 0.0f,
			u, v));
	}

	for (uint32 i = 0; i < sliceCount; ++i)
	{
		// Winding so normal points downward
		meshData.Indices32.push_back(baseCenterIndex);
		meshData.Indices32.push_back(baseStartIndex + i + 1);
		meshData.Indices32.push_back(baseStartIndex + i);
	}

	return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreatePyramid(float width, float height, float depth)
{
	MeshData meshData;

	float w2 = 0.5f * width;
	float d2 = 0.5f * depth;

	XMFLOAT3 p0(-w2, 0.0f, -d2);
	XMFLOAT3 p1(+w2, 0.0f, -d2);
	XMFLOAT3 p2(+w2, 0.0f, +d2);
	XMFLOAT3 p3(-w2, 0.0f, +d2);
	XMFLOAT3 apex(0.0f, height, 0.0f);

	auto AddTriFlat = [&](const XMFLOAT3& a, const XMFLOAT3& b, const XMFLOAT3& c)
		{
			XMVECTOR A = XMLoadFloat3(&a);
			XMVECTOR B = XMLoadFloat3(&b);
			XMVECTOR C = XMLoadFloat3(&c);

			XMVECTOR n = XMVector3Normalize(XMVector3Cross(B - A, C - A));
			XMFLOAT3 nf; XMStoreFloat3(&nf, n);

			// Tangent can be anything consistent for now
			XMFLOAT3 t(1.0f, 0.0f, 0.0f);

			uint32 start = (uint32)meshData.Vertices.size();
			meshData.Vertices.push_back(Vertex(a.x, a.y, a.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 0.0f, 1.0f));
			meshData.Vertices.push_back(Vertex(b.x, b.y, b.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 0.5f, 0.0f));
			meshData.Vertices.push_back(Vertex(c.x, c.y, c.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 1.0f, 1.0f));

			meshData.Indices32.push_back(start + 0);
			meshData.Indices32.push_back(start + 1);
			meshData.Indices32.push_back(start + 2);
		};

	// 4 side faces
	AddTriFlat(p0, p1, apex);
	AddTriFlat(p1, p2, apex);
	AddTriFlat(p2, p3, apex);
	AddTriFlat(p3, p0, apex);

	// Base (two triangles, normal down)
	uint32 base = (uint32)meshData.Vertices.size();
	meshData.Vertices.push_back(Vertex(p0.x, p0.y, p0.z, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f));
	meshData.Vertices.push_back(Vertex(p1.x, p1.y, p1.z, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f));
	meshData.Vertices.push_back(Vertex(p2.x, p2.y, p2.z, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f));
	meshData.Vertices.push_back(Vertex(p3.x, p3.y, p3.z, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f));

	meshData.Indices32.push_back(base + 0);
	meshData.Indices32.push_back(base + 2);
	meshData.Indices32.push_back(base + 1);

	meshData.Indices32.push_back(base + 0);
	meshData.Indices32.push_back(base + 3);
	meshData.Indices32.push_back(base + 2);

	return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateWedge(float width, float height, float depth)
{
	MeshData meshData;

	float w2 = 0.5f * width;
	float d2 = 0.5f * depth;

	// Bottom rectangle
	XMFLOAT3 b0(-w2, 0.0f, -d2);
	XMFLOAT3 b1(+w2, 0.0f, -d2);
	XMFLOAT3 b2(+w2, 0.0f, +d2);
	XMFLOAT3 b3(-w2, 0.0f, +d2);

	// Top is sloped: front is tall, back is shorter
	XMFLOAT3 t0(-w2, height, -d2);
	XMFLOAT3 t1(+w2, height, -d2);
	XMFLOAT3 t2(+w2, 0.5f * height, +d2);
	XMFLOAT3 t3(-w2, 0.5f * height, +d2);

	auto AddQuadFlat = [&](const XMFLOAT3& a, const XMFLOAT3& b, const XMFLOAT3& c, const XMFLOAT3& d)
		{
			XMVECTOR A = XMLoadFloat3(&a);
			XMVECTOR B = XMLoadFloat3(&b);
			XMVECTOR C = XMLoadFloat3(&c);

			XMVECTOR n = XMVector3Normalize(XMVector3Cross(B - A, C - A));
			XMFLOAT3 nf; XMStoreFloat3(&nf, n);

			XMFLOAT3 t(1, 0, 0);

			uint32 start = (uint32)meshData.Vertices.size();
			meshData.Vertices.push_back(Vertex(a.x, a.y, a.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 0, 1));
			meshData.Vertices.push_back(Vertex(b.x, b.y, b.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 1, 1));
			meshData.Vertices.push_back(Vertex(c.x, c.y, c.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 1, 0));
			meshData.Vertices.push_back(Vertex(d.x, d.y, d.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 0, 0));

			meshData.Indices32.push_back(start + 0); meshData.Indices32.push_back(start + 1); meshData.Indices32.push_back(start + 2);
			meshData.Indices32.push_back(start + 0); meshData.Indices32.push_back(start + 2); meshData.Indices32.push_back(start + 3);
		};

	// Faces
	AddQuadFlat(b0, b1, b2, b3);     // bottom
	AddQuadFlat(t0, t3, t2, t1);     // top (sloped)
	AddQuadFlat(b0, b1, t1, t0);     // front
	AddQuadFlat(b3, b2, t2, t3);     // back
	AddQuadFlat(b0, t0, t3, b3);     // left
	AddQuadFlat(b1, b2, t2, t1);     // right

	return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateDiamond(float width, float height, float depth)
{
	MeshData meshData;

	float w2 = 0.5f * width;
	float d2 = 0.5f * depth;
	float h2 = 0.5f * height;

	XMFLOAT3 top(0.0f, +h2, 0.0f);
	XMFLOAT3 bot(0.0f, -h2, 0.0f);

	XMFLOAT3 p0(-w2, 0.0f, -d2);
	XMFLOAT3 p1(+w2, 0.0f, -d2);
	XMFLOAT3 p2(+w2, 0.0f, +d2);
	XMFLOAT3 p3(-w2, 0.0f, +d2);

	auto AddTriFlat = [&](const XMFLOAT3& a, const XMFLOAT3& b, const XMFLOAT3& c)
		{
			XMVECTOR A = XMLoadFloat3(&a);
			XMVECTOR B = XMLoadFloat3(&b);
			XMVECTOR C = XMLoadFloat3(&c);

			XMVECTOR n = XMVector3Normalize(XMVector3Cross(B - A, C - A));
			XMFLOAT3 nf; XMStoreFloat3(&nf, n);

			XMFLOAT3 t(1, 0, 0);

			uint32 start = (uint32)meshData.Vertices.size();
			meshData.Vertices.push_back(Vertex(a.x, a.y, a.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 0, 1));
			meshData.Vertices.push_back(Vertex(b.x, b.y, b.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 0.5f, 0));
			meshData.Vertices.push_back(Vertex(c.x, c.y, c.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 1, 1));

			meshData.Indices32.push_back(start + 0);
			meshData.Indices32.push_back(start + 1);
			meshData.Indices32.push_back(start + 2);
		};

	// Top half
	AddTriFlat(p0, p1, top);
	AddTriFlat(p1, p2, top);
	AddTriFlat(p2, p3, top);
	AddTriFlat(p3, p0, top);

	// Bottom half (flip winding)
	AddTriFlat(p1, p0, bot);
	AddTriFlat(p2, p1, bot);
	AddTriFlat(p3, p2, bot);
	AddTriFlat(p0, p3, bot);

	return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateTriangularPrism(float width, float height, float depth)
{
	MeshData meshData;

	float w2 = 0.5f * width;
	float d2 = 0.5f * depth;

	// Triangle (front z=-d2) and (back z=+d2)
	XMFLOAT3 f0(-w2, 0.0f, -d2);
	XMFLOAT3 f1(+w2, 0.0f, -d2);
	XMFLOAT3 f2(0.0f, height, -d2);

	XMFLOAT3 b0(-w2, 0.0f, +d2);
	XMFLOAT3 b1(+w2, 0.0f, +d2);
	XMFLOAT3 b2(0.0f, height, +d2);

	auto AddTriFlat = [&](const XMFLOAT3& a, const XMFLOAT3& b, const XMFLOAT3& c)
		{
			XMVECTOR A = XMLoadFloat3(&a);
			XMVECTOR B = XMLoadFloat3(&b);
			XMVECTOR C = XMLoadFloat3(&c);

			XMVECTOR n = XMVector3Normalize(XMVector3Cross(B - A, C - A));
			XMFLOAT3 nf; XMStoreFloat3(&nf, n);

			XMFLOAT3 t(1, 0, 0);

			uint32 start = (uint32)meshData.Vertices.size();
			meshData.Vertices.push_back(Vertex(a.x, a.y, a.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 0, 1));
			meshData.Vertices.push_back(Vertex(b.x, b.y, b.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 1, 1));
			meshData.Vertices.push_back(Vertex(c.x, c.y, c.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 0.5f, 0));

			meshData.Indices32.push_back(start + 0);
			meshData.Indices32.push_back(start + 1);
			meshData.Indices32.push_back(start + 2);
		};

	auto AddQuadFlat = [&](const XMFLOAT3& a, const XMFLOAT3& b, const XMFLOAT3& c, const XMFLOAT3& d)
		{
			XMVECTOR A = XMLoadFloat3(&a);
			XMVECTOR B = XMLoadFloat3(&b);
			XMVECTOR C = XMLoadFloat3(&c);

			XMVECTOR n = XMVector3Normalize(XMVector3Cross(B - A, C - A));
			XMFLOAT3 nf; XMStoreFloat3(&nf, n);

			XMFLOAT3 t(1, 0, 0);

			uint32 start = (uint32)meshData.Vertices.size();
			meshData.Vertices.push_back(Vertex(a.x, a.y, a.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 0, 1));
			meshData.Vertices.push_back(Vertex(b.x, b.y, b.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 1, 1));
			meshData.Vertices.push_back(Vertex(c.x, c.y, c.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 1, 0));
			meshData.Vertices.push_back(Vertex(d.x, d.y, d.z, nf.x, nf.y, nf.z, t.x, t.y, t.z, 0, 0));

			meshData.Indices32.push_back(start + 0); meshData.Indices32.push_back(start + 1); meshData.Indices32.push_back(start + 2);
			meshData.Indices32.push_back(start + 0); meshData.Indices32.push_back(start + 2); meshData.Indices32.push_back(start + 3);
		};

	// Front and back triangles
	AddTriFlat(f0, f2, f1);
	AddTriFlat(b0, b1, b2);

	// 3 rectangular sides
	AddQuadFlat(f0, f1, b1, b0); // bottom
	AddQuadFlat(f0, b0, b2, f2); // left
	AddQuadFlat(f1, f2, b2, b1); // right

	return meshData;
}

GeometryGenerator::MeshData GeometryGenerator::CreateTorus(float majorRadius, float minorRadius, uint32 sliceCount, uint32 stackCount)
{
	MeshData meshData;

	// i: around major circle (theta), j: around tube (phi)
	for (uint32 i = 0; i <= sliceCount; ++i)
	{
		float u = (float)i / (float)sliceCount;
		float theta = u * XM_2PI;

		float cosT = cosf(theta);
		float sinT = sinf(theta);

		// Tube center at this theta
		XMFLOAT3 center(majorRadius * cosT, 0.0f, majorRadius * sinT);

		for (uint32 j = 0; j <= stackCount; ++j)
		{
			float v = (float)j / (float)stackCount;
			float phi = v * XM_2PI;

			float cosP = cosf(phi);
			float sinP = sinf(phi);

			float xLocal = minorRadius * cosP;
			float yLocal = minorRadius * sinP;

			float x = (majorRadius + xLocal) * cosT;
			float y = yLocal;
			float z = (majorRadius + xLocal) * sinT;

			// Normal = from center to surface point
			XMFLOAT3 p(x, y, z);
			XMVECTOR P = XMLoadFloat3(&p);
			XMVECTOR C = XMLoadFloat3(&center);
			XMVECTOR n = XMVector3Normalize(P - C);

			XMFLOAT3 nf; XMStoreFloat3(&nf, n);

			// Tangent around theta
			XMFLOAT3 tangent(-sinT, 0.0f, cosT);

			meshData.Vertices.push_back(Vertex(
				x, y, z,
				nf.x, nf.y, nf.z,
				tangent.x, tangent.y, tangent.z,
				u, 1.0f - v));
		}
	}

	uint32 ring = stackCount + 1;
	for (uint32 i = 0; i < sliceCount; ++i)
	{
		for (uint32 j = 0; j < stackCount; ++j)
		{
			uint32 a = i * ring + j;
			uint32 b = (i + 1) * ring + j;
			uint32 c = (i + 1) * ring + (j + 1);
			uint32 d = i * ring + (j + 1);

			meshData.Indices32.push_back(a);
			meshData.Indices32.push_back(b);
			meshData.Indices32.push_back(c);

			meshData.Indices32.push_back(a);
			meshData.Indices32.push_back(c);
			meshData.Indices32.push_back(d);
		}
	}

	return meshData;
}
