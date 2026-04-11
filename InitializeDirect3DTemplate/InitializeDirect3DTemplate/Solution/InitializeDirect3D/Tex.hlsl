Texture2D gDiffuseMap : register(t0);

SamplerState gsamPointWrap : register(s0);
SamplerState gsamPointClamp : register(s1);
SamplerState gsamLinearWrap : register(s2);
SamplerState gsamLinearClamp : register(s3);
SamplerState gsamAnisotropicWrap : register(s4);
SamplerState gsamAnisotropicClamp : register(s5);

// ==========================
// PASS CONSTANTS (ONLY ONE)
// ==========================
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

    // Directional/Ambient lighting
    float4 gAmbientLight;

    float3 gLightDir;
    float pad1;

    float3 gLightStrength;
    float pad2;

    // Point light
    float3 gPointLightPosition;
    float gPointLightFalloffStart;

    float3 gPointLightStrength;
    float gPointLightFalloffEnd;
};

// ==========================
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

// ==========================
struct VertexIn
{
    float3 PosL : POSITION;
    float3 NormalL : NORMAL;
    float2 TexC : TEXCOORD;
};

struct VertexOut
{
    float4 PosH : SV_POSITION;
    float3 PosW : POSITION1;
    float3 Normal : NORMAL;
    float2 TexC : TEXCOORD;
};

VertexOut VS(VertexIn vin)
{
    VertexOut vout;

    float4 posW = mul(float4(vin.PosL, 1.0f), gWorld);
    vout.PosW = posW.xyz;
    vout.PosH = mul(posW, gViewProj);

    vout.Normal = vin.NormalL;

    float4 texC = mul(float4(vin.TexC, 0.0f, 1.0f), gTexTransform);
    texC = mul(texC, gMatTransform);
    vout.TexC = texC.xy;

    return vout;
}

float4 PS(VertexOut pin) : SV_Target
{
    float3 normal = normalize(pin.Normal);

    // Directional light
    float3 dirLightDir = normalize(-gLightDir);
    float dirDiffuse = saturate(dot(normal, dirLightDir));

    float3 lighting = gAmbientLight.rgb +
                      dirDiffuse * gLightStrength;

    // Point light
    float3 toLight = gPointLightPosition - pin.PosW;
    float distanceToLight = length(toLight);
    float3 pointLightDir = normalize(toLight);

    float pointDiffuse = saturate(dot(normal, pointLightDir));

    float attenuation = saturate(
        (gPointLightFalloffEnd - distanceToLight) /
        (gPointLightFalloffEnd - gPointLightFalloffStart)
    );

    lighting += pointDiffuse * gPointLightStrength * attenuation;

    float4 texColor = gDiffuseMap.Sample(gsamAnisotropicWrap, pin.TexC);
    float4 diffuseAlbedo = texColor * gDiffuseAlbedo;

    float3 finalColor = diffuseAlbedo.rgb * lighting;

    return float4(finalColor, diffuseAlbedo.a);
}