#include "hg_sdf.glsl"

#define PRECISION 1000.
#define MAX_ITERATIONS 200
#define MAX_DIST 10.
#define EPSILON (1. / PRECISION + 0.0001)
#define time (iTime * (135. / 60.) / 4.)
const vec3 purple = normalize(vec3(0.298, 0.176, 0.459));

// https://thebookofshaders.com/11/
float rand (in vec2 st) {
    return fract(sin(dot(st.xy, vec2(12.9898, 78.233))) * 43758.5453123);
}

// 2D Noise based on Morgan McGuire @morgan3d
// https://www.shadertoy.com/view/4dS3Wd
float noise (in vec2 st) {
    vec2 i = floor(st);
    vec2 f = fract(st);
    float a = rand(i);
    float b = rand(i + vec2(1.0, 0.0));
    float c = rand(i + vec2(0.0, 1.0));
    float d = rand(i + vec2(1.0, 1.0));
    vec2 u = smoothstep(0.,1.,f);
    return mix(a, b, u.x) + (c - a) * u.y * (1.0 - u.x) + (d - b) * u.x * u.y;
}

void pp(inout vec3 c, vec2 uv)
{
    c = clamp(c, 0.02, 1.);
    c -= 0.02;
    c *= 1.1;
    c = sqrt(c);
    c = c*c*(2.5-1.45*c*c); // contrast
    c = pow(c, vec3(1.0,0.96,1.0)); // soft green
    c *= vec3(1.08,0.99,0.99); // tint red
    c.b = (c.b+0.05)/1.05; // bias blue
    // c = mix(c, c.ggg, 0.12); // desaturate
    c = 0.06 + (c * 0.9);
    
    c -= smoothstep(0.55, 1.3, abs(uv.x)) * 0.2; // vignette x
    c -= smoothstep(0.17, 0.7, abs(uv.y)) * 0.2; // vignette y
    c -= step(0.515, abs(uv.x)); // letterbox x
    c -= step(0.35, abs(uv.y)); // letterbox y
    c += rand(uv.xy + iTime) * 0.05; // noise
}

float displacement(vec3 p, float dist)
{
	return sin(dist*p.x)*sin(dist*p.y)*sin(dist*p.z);
}

float glass(float x, float sharpness) {
    return (1. / max(0., x)) / (sharpness * 100.);
}

float lightness(float df, float density, float sharpness)
{
    return min(1., glass(df, sharpness)) * density * 10.;
}

void add(inout float d, inout vec3 c, float d2, vec3 color, float density, float sharpness)
{
    d = min(d, d2);
    c += lightness(d2, density, sharpness) * color;
}

void sub(inout float d, inout vec3 c, float d2, vec3 color, float density, float sharpness)
{
    d = max(-d,d2);
    c -= lightness(d2, density, sharpness) * color;
}

vec4 scene(vec3 p)
{
    float n = time * 1. + 4.;
    float d = 1. / 0.;
    vec3 c = vec3(0.);

    p.z -= 0.5;

    pR(p.yz, 2.1);
    pR(p.xy, PI/2. + sin(n) * 1.5);
    add(d, c, max(-fDodecahedron(p, 0.15), fDodecahedron(p, 0.15)), purple, 3., 20.);
    pR(p.yx,  n);
    p.z += tan(n * 2.) * 0.1;
    add(d, c, max(-fBox(p, vec3(0.2, 0.2, 0.2)), fBox(p, vec3(0.2, 0.2, 0.04))), purple.grb, 2., 0.4);

    return vec4(c, d);
}

vec3 normalAt(in vec3 p)
{
    float c = scene(p).w;
	vec2 nOfs = vec2(0.001, 0.0);
	return normalize(vec3(scene(p + nOfs.xyy).w, scene(p + nOfs.yxy).w, scene(p + nOfs.yyx).w) - c);
}

vec3 march(vec3 pos, vec3 rayDir)
{
    vec3 c = vec3(0.);
    int reflections = 0;
    
    for(int i = 0; i < MAX_ITERATIONS; i++)
    {
        vec4 s = scene(pos);
        float d = max(s.w, (1. + (0.1 * pos.z)) / PRECISION);
        
        pos += rayDir * d;
        c += s.rgb / (20. + PRECISION);
        
        if(d < EPSILON && reflections < 1 )
        {
            reflections++;

            vec3 normal = normalAt(pos - rayDir * EPSILON * 2.);
            vec3 rayDir = reflect(rayDir, normal);
            vec3 pos = pos;

            for(int i = 0; i < MAX_ITERATIONS; i++)
            {
                vec4 s = scene(pos);
                float d = max(s.w, (1. + (0.1 * pos.z)) / PRECISION);
                
                pos += rayDir * d;
                c += (s.rgb / PRECISION) * 0.3;

                if(d > MAX_DIST)
                    break;
            }
        }
        else if(d > MAX_DIST)
            break;
        else if(i == MAX_ITERATIONS - 1)
            c = vec3(1., 0., 0.);
    }

    return c;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (fragCoord.xy - (iResolution.xy * 0.5)) / iResolution.yy;
    vec3 origin = vec3(0., 0., -1.);
    vec3 cameraTarget = origin + vec3(0., 0., 1.);
    vec3 upDirection = vec3(0., 1.0, 0.);
    vec3 cameraDir = normalize(cameraTarget - origin);
    vec3 cameraRight = normalize(cross(upDirection, origin));
	vec3 cameraUp = cross(cameraDir, cameraRight);
    vec3 rayDir = normalize(cameraRight * uv.x + cameraUp * uv.y + cameraDir);

    vec3 c = march(origin, rayDir);
    pp(c, uv);
   
    fragColor = vec4(c, 1.);
}