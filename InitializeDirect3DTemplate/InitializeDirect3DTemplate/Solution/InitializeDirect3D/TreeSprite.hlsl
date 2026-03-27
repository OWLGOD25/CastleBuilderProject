Texture2D gDiffuseMap : register(t0);

SamplerState gsamPointWrap : register(s0);
SamplerState gsamPointClamp : register(s1);
SamplerState gsamLinearWrap : register(s2);
SamplerState gsamLinearClamp : register(s3);
SamplerState gsamAnisotropicWrap : register(s4);
SamplerState gsamAnisotropicClamp : register(s5);

cbuffer cbPass : register(b1)
{
    float4x4 gView;
    float4x4 gInvView;
    float4x4 gProj;
    float4x4 gInvProj;
    float4x4 gViewProj;
    float4x4 gInvViewProj;

    float3 gEyePosW;
    float cbPerObjectPad1;

    float2 gRenderTargetSize;
    float2 gInvRenderTargetSize;

    float gNearZ;
    float gFarZ;
    float gTotalTime;
    float gDeltaTime;

    float4 gAmbientLight;

    float3 gLightDir;
    float pad1;

    float3 gLightStrength;
    float pad2;
};

cbuffer cbPerObject : register(b0)
{
    float4x4 gWorld;
    float4x4 gTexTransform;
};

cbuffer cbMaterial : register(b2)
{
    float4 gDiffuseAlbedo;
    float3 gFresnelR0;
    float gRoughness;
    float4x4 gMatTransform;
};

struct VertexIn
{
    float3 PosW : POSITION;
    float2 SizeW : SIZE;
};

struct VertexOut
{
    float3 CenterW : POSITION;
    float2 SizeW : SIZE;
};

VertexOut VS(VertexIn vin)
{
    VertexOut vout;
    vout.CenterW = vin.PosW;
    vout.SizeW = vin.SizeW;
    return vout;
}

struct GeoOut
{
    float4 PosH : SV_POSITION;
    float2 TexC : TEXCOORD;
};

[maxvertexcount(4)]
void GS(point VertexOut gin[1], inout TriangleStream<GeoOut> triStream)
{
    float3 center = gin[0].CenterW;
    float2 size = gin[0].SizeW;

    float3 up = float3(0.0f, 1.0f, 0.0f);
    float3 look = gEyePosW - center;
    look.y = 0.0f;
    look = normalize(look);

    float3 right = normalize(cross(up, look));

    float halfWidth = 0.5f * size.x;
    float halfHeight = 0.5f * size.y;

    float3 v[4];
    v[0] = center + halfWidth * right - halfHeight * up;
    v[1] = center + halfWidth * right + halfHeight * up;
    v[2] = center - halfWidth * right - halfHeight * up;
    v[3] = center - halfWidth * right + halfHeight * up;

    float2 tex[4] =
    {
        float2(1.0f, 1.0f),
        float2(1.0f, 0.0f),
        float2(0.0f, 1.0f),
        float2(0.0f, 0.0f)
    };

    GeoOut gout;

    gout.PosH = mul(float4(v[0], 1.0f), gViewProj);
    gout.TexC = tex[0];
    triStream.Append(gout);

    gout.PosH = mul(float4(v[1], 1.0f), gViewProj);
    gout.TexC = tex[1];
    triStream.Append(gout);

    gout.PosH = mul(float4(v[2], 1.0f), gViewProj);
    gout.TexC = tex[2];
    triStream.Append(gout);

    gout.PosH = mul(float4(v[3], 1.0f), gViewProj);
    gout.TexC = tex[3];
    triStream.Append(gout);
}

float4 PS(GeoOut pin) : SV_Target
{
    float4 c = gDiffuseMap.Sample(gsamLinearClamp, pin.TexC);
    clip(c.a - 0.1f);
    return c;
}