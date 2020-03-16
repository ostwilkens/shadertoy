#define PI 3.14159265358979323846
const vec3 purple = normalize(vec3(0.298, 0.176, 0.459));

void pR(inout vec2 p, float a) {
    // thx to hg (http://mercury.sexy/hg_sdf)
    p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

void pR45(inout vec2 p) {
    // thx to hg (http://mercury.sexy/hg_sdf)
	p = (p + vec2(p.y, -p.x))*sqrt(0.5);
}

// float lightField(float distanceField) {
//     return min(1., (1. / distanceField) * 0.005) * 0.2;
// }


float lightField(float distanceField) {
    return (1. / distanceField) * 0.0001;
}

float cube(vec3 p, float r) {
    return length(p - clamp(p, -r, r));
}

vec4 scene(vec3 p)
{
    p.z -= 0.1;

    pR45(p.xy);
    pR45(p.yz);
    pR45(p.zx);

    float innerCubeDist = cube(p, 0.048);
    float outerCubeDist = cube(p, 0.05);
    float d = outerCubeDist;

    float innerCubeLight = lightField(innerCubeDist);
    float outerCubeLight = lightField(outerCubeDist);
    float lf = outerCubeLight - innerCubeLight;

    vec3 l = vec3(0.0);
    l += lf * purple;

    return vec4(l, d);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (fragCoord.xy - (iResolution.xy * 0.5)) / iResolution.yy;

    vec3 cameraOrigin = vec3(0., 0., -0.5);
    vec3 cameraTarget = cameraOrigin + vec3(0., 0., 1.);
    vec3 upDirection = vec3(0., 1.0, 0.);
    vec3 cameraDir = normalize(cameraTarget - cameraOrigin);
    vec3 cameraRight = normalize(cross(upDirection, cameraOrigin));
	vec3 cameraUp = cross(cameraDir, cameraRight);
    vec3 rayDir = normalize(cameraRight * uv.x + cameraUp * uv.y + cameraDir);

    vec3 c = vec3(0.0);
    vec3 pos = cameraOrigin;

    for(int i = 0; i < 10; i++)
    {
        vec4 ld = scene(pos);
        vec3 l = ld.rgb;
        float dist = ld.w;

        pos += rayDir * dist;
        c += l;
    }
   
    fragColor = vec4(c, 1.);
}