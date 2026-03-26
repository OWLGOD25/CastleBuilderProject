cbuffer cbPass : register(b1)
{
	float3 gLightDir;
	float pad0;

	float3 gLightStrength;
	float pad1;

	float4 gAmbientLight;
};

struct VertexOut
{
	float4 PosH   : SV_POSITION;
	float4 Color  : COLOR;
	float3 Normal : NORMAL;
};

float4 PS(VertexOut pin) : SV_Target
{
	float3 normal = normalize(pin.Normal);

	float3 lightDir = normalize(-gLightDir);

	float diffuse = saturate(dot(normal, lightDir));

	float3 lighting = gAmbientLight.rgb +
					  diffuse * gLightStrength;

	float3 finalColor = pin.Color.rgb * lighting;

	return float4(finalColor, pin.Color.a);
}