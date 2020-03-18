#iChannel0 "file://noise.png"
#iChannel0::WrapMode "Repeat"
#include "hg_sdf.glsl"

#define PRECISION 1000.
#define MAX_ITERATIONS 250
#define MAX_DIST 10.
#define EPSILON (1. / PRECISION + 0.0001)
const vec3 purple = normalize(vec3(0.298, 0.176, 0.459));

float txnoise(in vec3 x)
{
    // The MIT License
    // Copyright © 2013 Inigo Quilez
    // Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY distIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);
    vec2 uv = (p.xy+vec2(37.0,17.0)*p.z) + f.xy;
    vec2 rg = textureLod( iChannel0, (uv+0.5)/256.0, 0.0).yx;
    return mix( rg.x, rg.y, f.z );
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
    return min(1., glass(df, sharpness)) * density * 7.;
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
    float n = iTime * 1.+ 4.;
    float d = 1. / 0.;
    vec3 c = vec3(0.);

    p.z -= 0.5;


    // p.x += 0.5;
    // pR(p.xz, sin(n) * 0.5);
    add(d, c, fBoxCheap(p + vec3(0., 0., -3.), vec3(1., 1., 0.01)), purple.gbr, 1., 10.);


    pR(p.yz, 2.);
    pR(p.xy, PI/2. + sin(n) * 1.5);
    // pR(p.yz, + n * 0.2);
    // pR(p.zx, + n);

    
    // add(d, c, max(-fSphere(p, 0.15), fSphere(p, 0.15)), purple, 2., 5.);
    add(d, c, max(-fDodecahedron(p, 0.15), fDodecahedron(p, 0.15)), purple, 5., 20.);


    // p.x += 0.4;

    // add(d, c, fBox(p, vec3(0.1, 0.1, 0.1)), vec3(0., 1., 0.), 1., 1.);


    pR(p.yx,  n);
    p.z += tan(n * 2.) * 0.1;
    add(d, c, max(-fBox(p, vec3(0.2, 0.2, 0.2)), fBox(p, vec3(0.2, 0.2, 0.04))), purple.grb, 3., 0.4);

    // p.x -= 0.15;
    // add(d, c, fBox(p, vec3(0.1)), purple, 9., 1.);
    // sub(d, c, fBox(p, vec3(0.098)), purple, 9., 6.);

    // p.x += 0.3;
    // add(d, c, fBox(p, vec3(0.01, 0.1, 0.1)), purple.gbr, 1., 0.1);


    return vec4(c, d);
}

vec3 normalAt(in vec3 p)
{
    vec2 e = vec2(-EPSILON, EPSILON);
    float t1 = scene(p + e.yxx).w; 
    float t2 = scene(p + e.xxy).w;
    float t3 = scene(p + e.xyx).w;
    float t4 = scene(p + e.yyy).w;

    return vec3(normalize(e.yxx*t1 + e.xxy*t2 + e.xyx*t3 + e.yyy*t4));
}

vec3 normalAt2(in vec3 p)
{
	vec3 normal;
	vec3 ep = vec3(EPSILON, 0, 0);

	normal.x = scene(p + ep.xyz).w - scene(p - ep.xyz).w;
	normal.y = scene(p + ep.yxz).w - scene(p - ep.yxz).w;
	normal.z = scene(p + ep.yzx).w - scene(p - ep.yzx).w;

	return normalize(normal);
}

vec3 normalAt3(in vec3 p)
{
    float c = scene(p).w;
	vec2 nOfs = vec2(0.001, 0.0);
	return normalize(vec3(scene(p + nOfs.xyy).w, scene(p + nOfs.yxy).w, scene(p + nOfs.yyx).w) - c);
}

vec3 marchReflection(vec3 pos, vec3 rayDir)
{
    vec3 c = vec3(0.);
    vec3 normal = normalAt3(pos - rayDir * EPSILON * 2.);
    rayDir = reflect(rayDir, normal);
    
    for(int i = 0; i < MAX_ITERATIONS; i++)
    {
        vec4 s = scene(pos);
        float d = max(s.w, (1. + (0.1 * pos.z)) / (PRECISION * 1.));
        
        pos += rayDir * d;
        c += s.rgb / PRECISION;
        
        if(d > MAX_DIST)
            break;
    }

    return c;
}

vec3 march(vec2 uv)
{
    vec3 origin = vec3(0., 0., -1.);
    vec3 cameraTarget = origin + vec3(0., 0., 1.);
    vec3 upDirection = vec3(0., 1.0, 0.);
    vec3 cameraDir = normalize(cameraTarget - origin);
    vec3 cameraRight = normalize(cross(upDirection, origin));
	vec3 cameraUp = cross(cameraDir, cameraRight);
    vec3 rayDir = normalize(cameraRight * uv.x + cameraUp * uv.y + cameraDir);

    vec3 c = vec3(0.);
    vec3 pos = origin;

    float reflected = 0.;

    for(int i = 0; i < MAX_ITERATIONS; i++)
    {
        vec4 s = scene(pos);
        float d = max(s.w, (1. + (0.1 * pos.z)) / PRECISION);

        pos += rayDir * d;
        c += s.rgb / (20. + PRECISION);

        if(reflected <= 0. && d < EPSILON)
        {
            reflected = 1.;

            c += marchReflection(pos, rayDir) * 0.5;

            // c.r = 1.;

            // break;

            // rayDir = reflect(rayDir, normalize(rayDir));
            // c = vec3(1.);
        }
        // else 
        // {
        //     reflected -= EPSILON * 100.;
        // }
        

        if(d > MAX_DIST)
            break;
        else if(i == MAX_ITERATIONS - 1)
            c = vec3(1., 0., 0.);
    }

    return c;
}

vec3 post(vec3 c, vec2 uv)
{
    // grading
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
    // c = clamp(c, 0., 1.);
    c = max(c, txnoise(vec3(uv.xy, 0.) * 700.) * 0.01); // texturize black
    c += txnoise(vec3(uv.xy, sin(cos(iTime) * 1000.)) * 1000.) * 0.03; // noise

    return c;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (fragCoord.xy - (iResolution.xy * 0.5)) / iResolution.yy;

    vec3 c = march(uv);
    c = post(c, uv);
   
    fragColor = vec4(c, 1.);
}