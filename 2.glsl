#iChannel0 "file://noise.png"
#iChannel0::WrapMode "Repeat"
#iChannel0::MagFilter "Linear"
#include "hg_sdf.glsl"

#define PRECISION 1000.
#define MAX_ITERATIONS 100
#define MAX_DIST 10.
#define EPSILON (1. / PRECISION + 0.0008)
#define time (iTime * (155. / 60.) / 1.)
const vec3 purple = normalize(vec3(0.298, 0.176, 0.459));

struct Result
{
    vec3 c;
    float d;
    float o;
};

float fTriPrism(vec3 p, vec2 h)
{
    const float k = sqrt(3.0);
    h.x *= 0.5*k;
    p.xy /= h.x;
    p.x = abs(p.x) - 1.0;
    p.y = p.y + 1.0/k;
    if( p.x+k*p.y>0.0 ) p.xy=vec2(p.x-k*p.y,-k*p.x-p.y)/2.0;
    p.x -= clamp( p.x, -2.0, 0.0 );
    float d1 = length(p.xy)*sign(-p.y)*h.x;
    float d2 = abs(p.z)-h.y;
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}

float txnoise(in vec3 x)
{
    // The MIT License
    // Copyright © 2013 Inigo Quilez
    // Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);
    vec2 uv = (p.xy+vec2(37.0,17.0)*p.z) + f.xy;
    vec2 rg = textureLod(iChannel0, (uv+0.5)/256.0, 0.0).yx;
    return mix( rg.x, rg.y, f.z );
}

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

float noise2 (in vec2 st) {
    return txnoise(vec3(st, 1.));
}

vec2 hash( vec2 x )  // replace this by something better
{
    const vec2 k = vec2( 0.3183099, 0.3678794 );
    x = x*k + k.yx;
    return -1.0 + 2.0*fract( 16.0 * k*fract( x.x*x.y*(x.x+x.y)) );
}

float noise3( in vec2 p )
{
    vec2 i = floor( p );
    vec2 f = fract( p );
	
	vec2 u = f*f*(3.0-2.0*f);

    return mix( mix( dot( hash( i + vec2(0.0,0.0) ), f - vec2(0.0,0.0) ), 
                     dot( hash( i + vec2(1.0,0.0) ), f - vec2(1.0,0.0) ), u.x),
                mix( dot( hash( i + vec2(0.0,1.0) ), f - vec2(0.0,1.0) ), 
                     dot( hash( i + vec2(1.0,1.0) ), f - vec2(1.0,1.0) ), u.x), u.y);
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
    // c -= step(0.515, abs(uv.x)); // letterbox x
    c -= step(0.35, abs(uv.y)); // letterbox y
    c += rand(uv.xy + iTime) * 0.02; // noise
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
    return min(1., glass(df, sharpness)) * density * (1. + sharpness / 10.) * (1. + (PRECISION / 500.)) * 0.3;
}

void add(inout float d, inout float o, inout vec3 c, float d2, vec3 color, float density, float sharpness, float opacity)
{
    d = min(d, d2);
    c += lightness(d2, density, sharpness) * color;
    o += lightness(d2, density, sharpness) * opacity;
}

void sub(inout float d, inout vec3 c, float d2, vec3 color, float density, float sharpness)
{
    d = max(-d,d2);
    c -= lightness(d2, density, sharpness) * color;
}

float shell(float d)
{
    return max(-d, d);
}

void scene3(vec3 p, float n, inout vec3 c, inout float d, inout float o)
{
    n *= 0.1;
    n = 0.5;

    int sidesA = int(mod(floor(n), 8.));
    int sidesB = int(mod(floor(n - 1.), 8.));

    float weight = smoothstep(0., 0.3, mod(n, 1.));

    float funky = weight + floor(n);
    pR(p.zy, -0.2);
    pR(p.zx, n * 0.5);
    pR(p.zx, funky * PI * 0.5);

    float shapeA = fGDF(p, 0.15, sidesA, sidesA + 5);
    float shapeB = fGDF(p, 0.15, sidesB, sidesB + 5);
    float shapeC = shapeA * weight + shapeB * (1. - weight);

    float expansion = sin(-n * 0.5) * 0.05 + 0.2;
    float shapeD = fGDF(p, 0.15 + expansion, sidesA, sidesA + 5);
    float shapeE = fGDF(p, 0.15 + expansion, sidesB, sidesB + 5);
    float shapeF = shapeD * weight + shapeE * (1. - weight);
    shapeF = shell(shapeF);
    float shapeG = fBox2Cheap(p.yx + tan(-n * PI + 1.9) * 0.05 - 0.05, vec2(0.02 + expansion * 0.05, 1.));
    float shapeH = max(shapeF, shapeG);
    
    add(d, o, c, shell(shapeC), purple, 6., 10., 0.5);
    add(d, o, c, shapeH, purple.grb, 9., 5., 1.0);
}

Result scene(vec3 p, float viewId)
{
    float n = time;
    vec3 c = vec3(0.);
    float d = 1. / 0.;
    float o = 0.0;

    n += + viewId * 11.0;

    scene3(p, n, c, d, o);

    return Result(c, d, o);
}

vec3 origin()
{
    return vec3(0., 0., -1.7);
}

vec3 target()
{
    return vec3(0., 0., 1.);
}

vec3 normalAt(in vec3 p, float viewId)
{
    float c = scene(p, viewId).d;
	vec2 nOfs = vec2(0.001, 0.0);
	return normalize(vec3(scene(p + nOfs.xyy, viewId).d, scene(p + nOfs.yxy, viewId).d, scene(p + nOfs.yyx, viewId).d) - c);
}

vec3 march(vec3 pos, vec3 rayDir, float viewId)
{
    vec3 c = vec3(0.);
    
    for(int i = 0; i < MAX_ITERATIONS; i++)
    {
        Result s = scene(pos, viewId);
        float d = max(s.d, (1. + (0.1 * pos.z)) / PRECISION);
        // float d = s.d + s.o;
        
        pos += rayDir * d;
        c += s.c / (20. + PRECISION);
        
    }

    return c;
}

vec3 image(vec2 uv, out vec4 fragColor, in vec2 fragCoord, float id)
{
    vec3 origin = origin();
    vec3 upDirection = vec3(0., 1.0, 0.);
    vec3 cameraDir = normalize(target() - origin);
    vec3 cameraRight = normalize(cross(upDirection, origin));
	vec3 cameraUp = cross(cameraDir, cameraRight);
    vec3 rayDir = normalize(cameraRight * uv.x + cameraUp * uv.y + cameraDir);

    vec3 c = march(origin, rayDir, id);

    return c;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (fragCoord.xy - (iResolution.xy * 0.5)) / iResolution.yy;
    vec3 c = vec3(0.0, 0.0, 0.0);

    c += image(uv, fragColor, fragCoord, 0.0);
    pp(c, uv);
    
    fragColor = vec4(c, 1.);
}