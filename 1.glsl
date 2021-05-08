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

float hash(float n) { return fract(sin(n) * 1e4); }

vec2 hash( vec2 x )  // replace this by something better
{
    const vec2 k = vec2( 0.3183099, 0.3678794 );
    x = x*k + k.yx;
    return -1.0 + 2.0*fract( 16.0 * k*fract( x.x*x.y*(x.x+x.y)) );
}

float noise4(vec3 x) {
	const vec3 step = vec3(110, 241, 171);
	vec3 i = floor(x);
	vec3 f = fract(x);
	float n = dot(i, step);
	vec3 u = f * f * (3.0 - 2.0 * f);
	return mix(mix(mix( hash(n + dot(step, vec3(0, 0, 0))), hash(n + dot(step, vec3(1, 0, 0))), u.x),
					mix( hash(n + dot(step, vec3(0, 1, 0))), hash(n + dot(step, vec3(1, 1, 0))), u.x), u.y),
				mix(mix( hash(n + dot(step, vec3(0, 0, 1))), hash(n + dot(step, vec3(1, 0, 1))), u.x),
					mix( hash(n + dot(step, vec3(0, 1, 1))), hash(n + dot(step, vec3(1, 1, 1))), u.x), u.y), u.z);
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

void add(inout Result r, float d2, vec3 color, float density, float sharpness)
{
    r.d = min(r.d, d2);
    r.c += lightness(d2, density, sharpness) * color;
}

void add(inout Result r, float d2, vec3 color, float density, float sharpness, float opacity)
{
    r.d = min(r.d, d2);
    r.c += lightness(d2, density, sharpness) * color;
    r.o += smoothstep(0.1, 0.0, d2) * opacity;
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

void scene2(vec3 p, float n, inout Result r)
{
    n *= 0.5;

    p.z += 1.0;

    // pR(p.xz, -0.4);
    // pR(p.xy, PI / 2.0);

    p.z -= 0.5;
    pR(p.yz, 2.1);
    pR(p.xy, PI/2. + sin(n) * 1.5);
    add(r, max(-fDodecahedron(p, 0.15), fDodecahedron(p, 0.15)), purple, 5., 10.);
    pR(p.yx,  n);
    p.z += tan(n * PI) * 0.1;
    add(r, max(-fBox(p, vec3(0.2, 0.2, 0.2)), fBox(p, vec3(0.2, 0.2, 0.04))), purple.grb, 5., 0.4);
}

void scene2b(vec3 p, float n, inout Result r)
{
    n *= 0.5;

    vec3 p1 = p;

    pR(p1.xz, -0.4);
    pR(p1.xy, PI / 2.0);

    p1.z -= 0.5;
    pR(p1.yz, 2.1);
    pR(p1.xy, PI/2. + sin(n) * 1.5);
    
    pR(p1.yx,  n);

    add(r, max(-fBox(p1, vec3(0.2, 0.2, 0.2)), fBox(p1, vec3(0.2, 0.2, 0.04))), purple.grb, 5., 0.4);

    p1.z += tan(n * PI) * 0.1;
    
    add(r, max(-fDodecahedron(p1, 0.15), fDodecahedron(p1, 0.15)), purple, 5., 10.);
}

void scene3(vec3 p, float n, inout Result r)
{
    // n *= 0.1;
    n *= 0.5;
    // n = 0.5;

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
    
    add(r, shell(shapeC), purple, 5., 10.);
    add(r, shapeH, purple.grb * 0.9, 10.0, 0.5);
}

float s4exp(float off, float n)
{
    return smoothstep(off, 0., n) / n * 0.1 + 
        smoothstep(off, 0., mod(n, 6.)) * smoothstep(0., 6., mod(n, 6.)) * 0.1 + 
        0.2 / (10. + mod(n, 1.) * 40.);
}

void pumpcube(vec3 p, inout Result r, float n, float s) 
{
    add(r, fBoxCheap(p + vec3(0., 0., s + s4exp(1., n)), vec3(s, s, 0.)), purple, 4., 20.);
    // p.y += tan(n) * 0.02;
    // pR(p.yz, n);
    add(r, fBoxCheap(p - vec3(0., 0., s + s4exp(2., n)), vec3(s, s, 0.)), purple, 4., 20.);
    // p.z += tan(n) * 0.02;
    // pR(p.zx, n);
    add(r, fBoxCheap(p + vec3(0., s + s4exp(3., n), 0.), vec3(s, 0., s)), purple, 4., 20.);
    // p.x += tan(n) * 0.02;
    // pR(p.xz, n);
    add(r, fBoxCheap(p - vec3(0., s + s4exp(4., n), 0.), vec3(s, 0., s)), purple, 4., 20.);
    // p.y += tan(n) * 0.02;
    // pR(p.yx, n);
    add(r, fBoxCheap(p + vec3(s + s4exp(5., n), 0., 0.), vec3(0., s, s)), purple, 4., 20.);
    // p.z += tan(n) * 0.02;
    // pR(p.zy, n);
    add(r, fBoxCheap(p - vec3(s + s4exp(5.5, n), 0., 0.), vec3(0., s, s)), purple, 4., 20.);
}

void scene4(vec3 p, float n, inout Result r)
{
    pR(p.zy, -0.2);
    pR(p.zx, n * 0.5);

    // p.x += tan(n) * 0.02;
    // pR(p.xy, n);

    pumpcube(p, r, n * 0.25, 0.1);
    // pR(p.zx, -n * 1.0);

    pumpcube(p, r, n * 0.25 + 0.25, 0.2 + sin(n * PI / 8.0) * 0.09);
}

void scene5(vec3 p, float n, inout Result r)
{
    float b = 1. - smoothstep(0., 0.3, mod(n, 1.));
    float bi = (1. - b) + floor(n);

    p.z += n * 0.5;

    add(r, fBox2Cheap(p.xy + vec2(0., -0.3), vec2(0.2, 0.)), purple.grb, 1., 10.);

    p.y += sin(p.z * PI) * 0.01;
    pR(p.xy, p.z + 1.5);

    add(r, fBox2Cheap(p.xy + vec2(0., 0.1), vec2(0., 0.003 + sin(p.z) * 0.004)), purple.grb, 10., 1.);

    p.y += noise(p.xz * 8. + n * 1.5) * 0.05;
    // p.y += txnoise(vec3(p.xz * 8. + n * 1.5, 0.)) * 0.05;
    p.x -= noise(p.xz * 2. + n * 0.5 + 10.) * 0.1;
    p.y += noise(p.xz * 0.5 - n * 0.25) * 0.15;
    add(r, fBox2Cheap(p.xy + vec2(0., 0.15), vec2(0.1, 0.)), purple, 5., 5.);
}

void scene6(vec3 p, float n, inout Result r)
{
    pR(p.yz, -0.2);
    vec3 p1 = p;

    pR(p.yz, 0.1 + p.z * 0.002);
    p.z += n * 1.;
    pR(p.xy, sin(p.z * 0.1) * 0.4);

    float f4 = noise(p.xz * 0.1) - 0.5;
    p.x += f4 * 4.;

    vec3 color = purple.grb * (1. - length(p1) * 0.008);

    add(r, shell(fBox2Cheap(p.xy + vec2(0., -100.05), vec2(0.05, 100.))), color, 2., 1.);
}

void scene7(vec3 p, float n, inout Result r) {
    p.z -= sin(n * 0.3) * 0.5 + 0.2;

    vec3 p1 = p;
    pR(p1.xy, n * 0.1);
    add(r, fBox2(p1.xy + vec2(-0.0001, -0.001), vec2(0., 0.)), purple.grb, 10., 10.);
    // add(r, fBox2Cheap(p.xy + vec2(0., -0.25), vec2(0.25, 0.)), purple.grb, 2., 10.);


    // pR(p.yz, n);
    // add(r, fCircle(p, 0.3) - 0.01, purple, 2., 10.);

    pR(p.xy, n * 0.2 + sin(n) * 0.1);
    float c1 = pModPolar(p.xy, 6.);

    p.x -= 0.3 + sin(c1) * 0.1;
    p.z += tan(n + c1) * 0.2;
    pR(p.zy, PI / 2.);
    pR(p.xz, n * 0.6 + sin(n + c1) * 0.1 + c1);
    // add(r, shell(fSphere(p, 0.08)), purple, 1., 2.);
    // pR(p.yz, sin(n * 0.5));
    // pR(p.yz, -0.5);
    pModPolar(p.xz, 3.);
    pR(p.xy, tan(n + c1 + 0.3));
    add(r, shell(fTriPrism(p, vec2(0.06, 0.007))), purple, 8., 5.);

}

void scene8(vec3 p, float n, inout Result r) {
    n *= 1.;
    n = mod(n, 122.);
    // n += 25.;
    // n = 55.;

    pR(p.yz, 0.2);
    p.y -= 0.3;

    vec3 p1 = p;
    vec3 p2 = p;

    p2.z -= 0.4;
    float height = -1. + n * 0.02;
    p2.y += height;
    float sun = (fSphere(p2, 0.7));

    p1.y += noise(p1.xz * 20. + n * 1.) * 0.012;
    p1.y += noise(p1.xz * 13. + n * 1.) * 0.015;
    p1.y += sin(p1.z * 3. + n * 0.2) * 0.006 * (0. + clamp(n * 0.25 - 2., 0., 24.));
    float ocean = fCylinder(p1, 1.05, 0.1);

    // ocean -= sun * p2.y * 0.1;
    ocean = fOpUnionRound(ocean, sun, max(0., 0.7 - n * 0.0065));
    ocean -= sun * 0.3;

    add(r, sun, purple.brg, 15., 0.3);
    add(r, ocean, purple.grb, 2., 5.);
}

float point(vec3 p, float r, float s) {
    return min(1., (1. / length(p)) * r) * s;
}

void scene9(vec3 p, float n, inout Result r)
{
    n *= 1.;
    n = mod(n, 10.);

    vec3 p1 = p;

    p1 += sin(p1.y * 10.0 + n) * 0.5;

    float sun = fSphere(p1, 0.2);

    add(r, sun, purple.brg, 25., 0.05);
}

void scene10(vec3 p, float n, inout Result r)
{
    // n *= 8.;
    n = mod(n, 30.);

    // n += 10.0;
    // n += 80.0;

    // n += 30.0;
    // n = 125.0;

    // float accel = max(0.0, n - 120.0);
    // float accel = smoothstep(117.0, 130.0, n);

    vec3 p1 = p;
    vec3 p2 = p;

    p1.x += sin(n * 0.25 - 1.0);
    p1.y += cos(n * 0.2);
    float sun = fSphere(p1, 0.2);
    // add(r, sun, purple.brg, 25., 0.05);
    // p1 += sin(p1.y * 10.0 + n) * 0.5;

    p2.y -= 0.2;
    pR(p2.yz, -0.5);
    pR(p2.xz, sin(n * 0.5) * 0.1);
    // p2.z += accel * 1.0 *;
    // p2.y -= accel * 1.0;
    // p2.y += 1.0 + accel;
    // pR(p2.yz, -0.5 - n * 0.041);
    pR(p2.xy, 0.9 + n * 0.02);
    pR(p2.yz, 0.5 + n * 0.0797);
    p.z += n;
    p2.xy *= 2.0;
    // p2.z += sin(p.y - 1.0) * 0.5;

    vec2 cell = pMod2(p2.xy, vec2(0.5));

    p2.z += (rand(cell.xy) - 0.5) * 0.5;

    float star = fSphere(p2, 0.005 + sin(cell.y) * 0.04);
    add(r, star, purple, 20.0, 0.1);
}

void scene11(vec3 p, float n, inout Result r)
{
    // n *= 8.;
    n = mod(n, 30.);

    // n += 10.0;
    // n += 80.0;

    // n += 30.0;
    // n = 125.0;

    vec3 p1 = p;
    vec3 p2 = p;

    // p1.x += sin(n * 0.25 - 1.0);
    // p1.y += cos(n * 0.2);
    float sun = fSphere(p1, 0.2);
    // add(r, sun, purple.brg, 25., 0.05);
    // p1 += sin(p1.y * 10.0 + n) * 0.5;

    // p2.y -= 0.8;
    // pR(p2.yz, -0.5);
    // pR(p2.xz, sin(n * 0.5) * 0.1);
    // p2.z += accel * 1.0 *;
    // p2.y -= accel * 1.0;
    // p2.y += 1.0 + accel;
    // pR(p2.yz, -0.5 - n * 0.041);
    // pR(p2.xy, 0.9 + n * 0.02);
    // p.z += n * 10.0;
    // pR(p2.yz, 0.5 + n * 0.0797);
    p2.xy += n;
    p2.xy *= 2.0;
    pR(p2.xz, sin(n * 0.2) * 0.1);
    pR(p2.zy, cos(n * 0.2) * 0.1);

    // p2.xy *= 1.25;

    // p2.xy *= sin(n * 0.25) + 2.0;
    p2.z += sin(p.y - 1.0) * 0.5;


    p2.z -= 2.0;

    vec2 cell = pMod2(p2.xy, vec2(0.5));
    p2.z /= 30.0;

    p2.z += (rand(cell.xy * 0.1) - 0.5) * 1.0;


    // float star = fSphere(p2, 0.005 + sin(cell.y) * 0.04);
    float star = fSphere(p2, 0.01);
    add(r, star, purple, 20.0, 0.1);
}


void scene12(vec3 p, float n, inout Result r)
{
    // n *= 8.;
    n = mod(n, 70.);

    // n = 53.0;


    vec3 p1 = p;
    vec3 p2 = p;

    // p1.x += sin(n * 0.25 - 1.0);
    // p1.y += cos(n * 0.2);
    float sun = fSphere(p1, 0.2);
    // add(r, sun, purple.brg, 25., 0.05);
    // p1 += sin(p1.y * 10.0 + n) * 0.5;

    p2.x -= 18.0;
    p2.z += pow(n, 1.7) * 0.1;

    // p2.y -= 0.8;
    pR(p2.yz, -0.5);
    // pR(p2.xz, sin(n * 0.5) * 0.1);

    
    pR(p2.yz, -smoothstep(0.0, 30.0, n) * 1.025); // panup++
    // pR(p2.yz, -1.0); // panup^
    

    
    pR(p2.xy, 0.9 + n * 0.01);
    // p.z += n * 10.0;
    // pR(p2.yz, 0.5 + n * 0.0797);
    // p2.xy += n;
    // p2.xy *= 2.0;
    // pR(p2.xz, sin(n * 0.2) * 0.1);
    pR(p2.zy, cos(n * 0.2) * 0.1);

    // p2.xy *= sin(n * 0.25) + 2.0;
    p2.z += sin(p.y - 1.0) * 0.5;


    p2.z -= 2.0;

    vec2 cell = pMod2(p2.xy, vec2(0.5));
    p2.z /= 3.0;

    // p2.z += (rand(cell.xy * 0.1) - 0.5) * 0.1;


    // float star = fSphere(p2, 0.005 + sin(cell.y) * 0.04);
    float star = fSphere(p2, 0.01);
    add(r, star, purple.grb, 22.0, 0.5);
}

float boxAbs(vec3 p, vec3 size) {
    return fBox(p + vec3(size.x, size.y, 0.0), size);
}

void scene13(vec3 p, float n, inout Result r)
{
    n = mod(n, 70.);
    // n = 0.0;

    // p.xy *= 1.5;

    // p.xy -= 0.5;
    // p.x -= 0.7;
    p.y += 0.1;
    p.y = -p.y;

    vec3 p1 = p;
    
    float glow = smoothstep(0.0, 0.9, 1.0 - 1. / n);
    float glow2 = smoothstep(0.6, 0.96, 1.0 - 1. / n);


    pR(p1.zx, noise(p1.xy * 2. + p1.z + n) * (1.0 - glow2));
    // pR(p1.xy, noise(p1.xz * 4. + n) * (1.0 - glow2));
    // pR(p1.xy, noise(p1.xz * 4. + n) - 0.5);


    pR(p1.xz, sin(n * 0.3) * 0.07);
    pR(p1.yz, cos(n * 0.42) * 0.05);
    p1.x -= 0.64;


    pR(p1.xz, 1.0 - glow2);
    p1.z -= 1.0 / (glow2 + 0.001);
    p1.z += 1.0;

    float maxdepth = 80.0;
    float depth = glow * maxdepth;
    p1.z += depth;
    p1.z -= maxdepth;


    pR(p1.xy, p1.z * 0.2);
    p1.y += sin(p1.z * PI * 0.22) * 0.1;
    // // add(r, fBox2Cheap(p.xy + vec2(0., 0.1), vec2(0., 0.003 + sin(p.z) * 0.004)), purple.grb, 10., 1.);
    // p1.y += noise(p1.zz * 1.0) * 0.1;
    // // p.y += txnoise(vec3(p.xz * 8. + n * 1.5, 0.)) * 0.05;
    p1.x -= noise(p1.xx * 2. + n * 0.2 + 10.) * 0.03;
    // p.y += noise(p.xz * 0.5 - n * 0.25) * 0.15;

    float a = boxAbs(p1, vec3(0.64, 0.17, maxdepth + .075 - depth + sin(n * 0.2) * 0.047 - n * 0.0018));
    a = max(-boxAbs(p1 + vec3(0.01, 0.01, 0.0), vec3(0.18, 0.07, 1000.0)), a);
    a = max(-boxAbs(p1 + vec3(0.01, 0.16, 0.0), vec3(0.18, 0.04, 1000.0)), a);
    a = max(-boxAbs(p1 + vec3(0.28, 0.01, 0.0), vec3(0.07, 0.3, 1000.0)), a);
    a = max(-boxAbs(p1 + vec3(0.54, 0.01, 0.0), vec3(0.07, 0.3, 1000.0)), a);
    a = max(-boxAbs(p1 + vec3(0.69, 0.01, 0.0), vec3(0.14, 0.07, 1000.0)), a);
    a = max(-boxAbs(p1 + vec3(0.69, 0.16, 0.0), vec3(0.14, 0.04, 1000.0)), a);
    a = max(-boxAbs(p1 + vec3(0.98, 0.01, 0.0), vec3(0.12, 0.04, 1000.0)), a);
    a = max(-boxAbs(p1 + vec3(0.98, 0.25, 0.0), vec3(0.07, 0.1, 1000.0)), a);
    a = max(-boxAbs(p1 + vec3(1.13, 0.16, 0.0), vec3(0.2, 0.2, 1000.0)), a);

    a -= 0.015;

    add(r, shell(a), purple, 100. - glow * 80.0, 1.5 - glow2);
}

void backdrop1(vec3 p, float n, inout Result r)
{
    n = mod(n, 70.);

    vec3 p1 = p;

    pR(p1.zy, 0.7);

    float a = fBox(p1, vec3(0.5, 0.3, 0.001));

    add(r, (a), purple.grb, 25., 0.1);
}

void backdrop2(vec3 p, float n, inout Result r)
{
    n = mod(n, 70.);

    vec3 p1 = p;

    p1.z += 1.0;

    // pR(p1.xz, PI / 2.0);



    // p1.y = abs(p1.y);

    // p1 = vec3(abs(p1.x), p1.y, p1.z);
    // p1.x -= 0.5;
    // pR(p1.zy, 0.0);
    // pR(p1.xz, sin(n));

    // float a = fBox(p1, vec3(0.5, 0.3, 0.001));
    // pR(p1.zx, PI / 4.0 * 1.0);
    // a = min(a, fBox(p1, vec3(0.5, 0.3, 0.001)));
    // pR(p1.zx, PI / 4.0 * 2.0);
    // a = min(a, fBox(p1, vec3(0.5, 0.3, 0.001)));
    // pR(p1.zx, PI / 4.0 * 3.0);
    // a = min(a, fBox(p1, vec3(0.5, 0.3, 0.001)));

    float a = fHexagonCircumcircle(p1, vec2(1.2, 2.0));
    a = shell(a);
    a = max(-fBoxCheap(p1 + vec3(0.0, 0.0, -0.5), vec3(0.5, 2.0, 0.5)), a);
    

    add(r, (a), purple.grb, 5., 0.5);
}

void backdrop3(vec3 p, float n, inout Result r)
{
    n = mod(n, 70.);

    vec3 p1 = p;

    pR(p1.yz, -1.15);

    p1.z -= 0.3;




    // p1.y = abs(p1.y);

    // p1 = vec3(abs(p1.x), p1.y, p1.z);
    // p1.x -= 0.5;
    // pR(p1.zy, 0.0);
    // pR(p1.xz, sin(n));

    // float a = fBox(p1, vec3(0.5, 0.3, 0.001));
    // pR(p1.zx, PI / 4.0 * 1.0);
    // a = min(a, fBox(p1, vec3(0.5, 0.3, 0.001)));
    // pR(p1.zx, PI / 4.0 * 2.0);
    // a = min(a, fBox(p1, vec3(0.5, 0.3, 0.001)));
    // pR(p1.zx, PI / 4.0 * 3.0);
    // a = min(a, fBox(p1, vec3(0.5, 0.3, 0.001)));

    float a = fBoxCheap(p1 , vec3(2.0, 1.5, 0.001));
    

    // add(r, (a), purple.grb, 5., 1.5, 0.0);
}

void scene14(vec3 p, float n, inout Result r)
{
    // pR(p.xz, n * 0.1);
    // pR(p.zy, n * 0.15);
    pR(p.xz, 0.1);
    pR(p.zy, -0.01);
    p *= 0.8;

    vec3 c = purple.brg;
    // c.g += sin(p.x * 10.0) * 10.0;
    // d += c.r;

    float dis1 = noise3(p.xy * p.z * 60.0 + n * 0.2);
    c += dis1 * 5.0 * c;
    p *= 1.0 - abs(dis1) * 0.1;

    // p.x /= 1.0 + smoothstep(0.1, 0.0, abs(p.y) * 0.5);
    p.y *= 1.5;
    p.xz *= 0.32 + clamp(0., 0.4, abs(p.y)) * 4.0;


    p.y *= 0.9;

    float d1 = fSphere(p, 0.48);

    add(r, d1, c, 5.0, 1.0, 1.0);
}


void scene15(vec3 p, float n, inout Result r)
{
    float turn = sin(n) * 0.4;
    // n -= 0.8;
    // n *= 0.5;


    // pR(p.xz, 1.1);
    n *= 2.0;
    


    vec3 p1 = p;
    
    p1.y += noise((p1.xz - n) * 0.2) * 0.2;
    // p1.y += noise((p1.xz - n) * 2.0) * 0.1;
    // p1.y -= noise((p1.xz - n) * 30.0) * 0.005;


    // pR(p1.yx, 0.05);
    pR(p1.xz, 0.7);

    // p.z -= 0.5;
    // pR(p.yz, 2.1);
    // pR(p.xy, PI/2. + sin(n) * 1.5);

    float d1 = fBox(p1, vec3(0.35, 0.2, 0.2));

    vec3 p2 = p1 + vec3(-0.4, 0.0, 0.0);
    pR(p2.xy, -0.15);
    d1 = max(-fBox(p2, vec3(0.3)), d1);
    
    vec3 p3 = p1 + vec3(0.36, -0.46, 0.0);
    pR(p3.xy, 0.17);
    d1 = max(-fBox(p3, vec3(0.3)), d1);
    
    vec3 p4 = p1 + vec3(0.23, 0.4, 0.0);
    pR(p4.xy, 0.6);
    d1 = max(-fBox(p4, vec3(0.3)), d1);

    vec3 p5 = p1 + vec3(0.5, 0.23, 0.0);
    pR(p5.xy, -0.17);
    d1 = max(-fBox(p5, vec3(0.3)), d1);

    vec3 p6 = p1 + vec3(0.0, 0., 0.59);
    pR(p6.yz, -0.1);
    d1 = max(-fBox(p6, vec3(0.4)), d1);
    
    vec3 p7 = p1 + vec3(0.0, 0., -0.59);
    pR(p7.yz, 0.1);
    d1 = max(-fBox(p7, vec3(0.4)), d1);

    add(r, shell(d1), vec3(1.0, 0.3, 0.2), 3., 1.0, 1.0);

    vec3 p8 = p1 + vec3(-0.22, -0.18, 0.0);
    float d2 = fBox(p8, vec3(0.6, 0.02, 0.22));

    vec3 p10 = p1 + vec3(-0.45, -0.1, 0.0);
    d2 = min(fBox(p10, vec3(0.3, 0.1, 0.2)), d2);
    
    add(r, shell(d2), vec3(1.0, 0.3, 0.05), 3., 1., 1.0);

    vec3 p9 = p1 + vec3(-0.22, -0.23, 0.0);
    float d3 = fBox(p9, vec3(0.45, 0.04, 0.2));
    
    vec3 p11 = p1 + vec3(-0.45, -0.03, 0.);
    d3 = min(d3, fBox(p11, vec3(0.3, 0.03, 0.2)));
    
    add(r, (d3), vec3(0.0), 1., 1., 1.0);

    
    vec3 p12 = p1 + vec3(0.3, -0.3, -0.23);
    pR(p12.xy, 0.05);
    pR(p12.xz, turn);
    float d4 = fBox(p12, vec3(0.2, 0.02, 0.05));

    
    vec3 p13 = p1 + vec3(0.3, -0.3, 0.23);
    pR(p13.xy, 0.05);
    pR(p13.xz, turn);
    d4 = min(d4, fBox(p13, vec3(0.2, 0.02, 0.05)));
    
    add(r, (d4), vec3(1.0, 0.0, 0.0), 2., 2., 1.0);

    
    vec3 p14 = p1 + vec3(0.15, 0.1, 0.);
    // pR(p14.xy, 0.05);
    // pR(p14.xz, turn);
    float d5 = fDodecahedron(p14, 0.05);
    
    add(r, d5, vec3(0.1, 0.65, 0.8) * 1., 10.0, 0.25, 1.0);
}


void scene16(vec3 p, float n, inout Result r)
{
    // n *= 0.2;
    n *= 0.5;

    pR(p.xz, n * 0.2);
    pR(p.zy, 0.7 + n * 0.1);

    // p1 += step(0.2, noise4(p * 10.0 + n)) * 0.1;

    float n1 = noise4(p * 1.7 + n * 0.25);

    p += n1 * 0.1;


    float d = fSphere(p, 0.4);

    vec3 p1 = p;

    p1 += n1;

    // p1.y *= sin(p.z);
    p1.x *= cos(p1.z * 3.) * 5.0;
    // p1.z *= sin(p1.x * 3.) * 5.0;
    // p1 /= tan(p1 * sin(n * 0.1)) + 1.0;

    pR(p1.xz, n * 0.2);
    p1.y = sin(p1.y);
    float i = 1.0 - step(0.1, noise4(p1 * 45.0));

    // d /= 1.0 + i * 10.0;
    // d -= i * 0.02;
    d -= i * 0.015;

    // vec3 p2 = p + 0.1;
    // pR(p2.xy, -n * 0.35);
    // p2.x = sin(p2.x);
    // p2.z = cos(p2.z);
    // i += 1.0 - step(0.03, noise4(p2 * 40.0));

    // vec3 p3 = p - 0.1;
    // pR(p3.yz, n * 0.17);
    // p3.x = sin(p2.x);
    // p3.z = cos(p2.z);
    // i += 1.0 - step(0.03, noise4(p3 * 40.0));


    // d = max(d, step(0.8, 1.0 - noise4(p * 10.0 + n) * 1.0) * 0.1);

    // d = max(d, step(1.0, noise4(p * 5.0) * 0.3));
    // d = max(d, noise4(0.4 * vec3(sin(p.x * 10.0), cos(p.y * 10.0), -abs(p.z) * 5.0)));

    add(r, (d), purple.grb, 20. * i, 20.0, 0.0);
}


void scene17(vec3 p, float n, inout Result r)
{
    // pR(p.xz, n * 0.2);


    pR(p.zy, 0.7 + n * 0.1);

    // p1 += step(0.2, noise4(p * 10.0 + n)) * 0.1;

    float d = fSphere(p, 0.5);

    vec3 p1 = p;

    // p1.y *= sin(p.z);
    p1.x *= cos(p1.z * 3.) * 5.0;
    p1.z *= sin(p1.x * 3.) * 5.0;
    // p1 *= tan(p1) + 1.0;

    pR(p1.xz, n * 0.2);
    // p1.y = sin(p1.y);
    // p1.z = cos(p1.z);
    float i = 1.0 - step(0.1, noise4(p1 * 40.0));

    // vec3 p2 = p + 0.1;
    // pR(p2.xy, -n * 0.35);
    // p2.x = sin(p2.x);
    // p2.z = cos(p2.z);
    // i += 1.0 - step(0.03, noise4(p2 * 40.0));

    // vec3 p3 = p - 0.1;
    // pR(p3.yz, n * 0.17);
    // p3.x = sin(p2.x);
    // p3.z = cos(p2.z);
    // i += 1.0 - step(0.03, noise4(p3 * 40.0));


    // d = max(d, step(0.8, 1.0 - noise4(p * 10.0 + n) * 1.0) * 0.1);

    // d = max(d, step(1.0, noise4(p * 5.0) * 0.3));
    // d = max(d, noise4(0.4 * vec3(sin(p.x * 10.0), cos(p.y * 10.0), -abs(p.z) * 5.0)));

    add(r, (d), purple.grb, 20. * i, 30.0, 0.0);
}


Result scene(vec3 p, float viewId)
{
    float n = time;
    // vec3 c = vec3(0.);
    // float d = 1. / 0.;
    // float o = 0.0;

    Result r = Result(vec3(0.), 1. / 0., 0.0);

    n += + viewId * 11.0;

    // scene2(p, n, r);
    // scene2b(p, n + 1.0, r);

    // backdrop2(p, n, r);
    // backdrop2(p, n, r);
    backdrop3(p, n, r);
    // scene3(p, n, r);

    // scene7(p, n, r);
    // scene8(p, n, r);

    // if(n > 8.0) {
        // scene13(p + vec3(0.0, 0.25, 0.0), n - 8.0, r);
    // }

    // scene13(p, n, r);
    // scene14(p, n, r);
    // scene15(p, n, r);
    scene16(p, n, r);
    
    return r;
}

vec3 origin1()
{
    return vec3(0., 0., -1.7 + sin(time) * 0.1);
}

vec3 origin()
{
    // return origin1();
    return vec3(0., 0., -1.7);
}

vec3 target1()
{
    // float rumbleX = noise(vec2(time * 10., 0.));
    // float rumbleY = noise(vec2(time * 10. + 100., 0.));
    // vec3 rumble = vec3(rumbleX, rumbleY, 0.);

    // float yoinkN = time * 2. + 0.7;
    // float yoinkX = sin(yoinkN * 70.);
    // vec3 yoink = vec3(yoinkX, 0., 0.) * max(0., (0.5 - fract(yoinkN)));
    
    return vec3(0., 0., 1.);
}

vec3 target()
{
    return target1();
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
    int reflections = 0;
    
    for(int i = 0; i < MAX_ITERATIONS; i++)
    {
        Result s = scene(pos, viewId);
        float d1 = max(s.d, (1. + (0.1 * pos.z)) / PRECISION);
        float d2 = s.d - 1.0 / PRECISION;
        float d = mix(d1, d2, s.o);
        
        pos += rayDir * d;
        c += s.c / (20. + PRECISION);
        
        vec3 normal = normalAt(pos + rayDir * EPSILON * 4., viewId);
        c *= 1.0 + (normal.x) * (normal.y) * (normal.z) * 0.05;

        // c /= 1.0 + d * 500.1;

        if(d < EPSILON && reflections < 1)
        {
            reflections++;

            vec3 normal = normalAt(pos - rayDir * EPSILON * 2., viewId);
            vec3 rayDir = reflect(rayDir, normal);
            vec3 pos = pos;

            for(int i = 0; i < MAX_ITERATIONS; i++)
            {
                Result s = scene(pos, viewId);
                float d = max(s.d, (1. + (0.1 * pos.z)) / PRECISION);
                
                pos += rayDir * d;
                c += (s.c / PRECISION) * 0.3;

                if(d > MAX_DIST)
                    break;
            }
        }
        else if(d > MAX_DIST)
            break;
        // else if(i == MAX_ITERATIONS - 1)
        //     c = vec3(1., 0., 0.);
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
    uv *= 1.5;

    // fragColor = vec4(1.0);
    vec3 c = vec3(0.0, 0.0, 0.0);

    c += image(uv, fragColor, fragCoord, 0.0);

    pp(c, uv);

    float h = 1.15;
    // if(uv.x < 0.0 && uv.y > 0.4) {
    //     c = image((uv * 2.0 + vec2(0.5, -h)), fragColor, fragCoord, 1.0);
    // } else if(uv.x > 0.0 && uv.y > 0.4) {
    //     c = image((uv * 2.0 + vec2(-0.5, -h)), fragColor, fragCoord, 2.0);
    // } else if(uv.x < 0.0 && uv.y < -0.4) {
    //     c = image((uv * 2.0 + vec2(0.5, h)), fragColor, fragCoord, 3.0);
    // } else if(uv.x > 0.0 && uv.y < -0.4) {
    //     c = image((uv * 2.0 + vec2(-0.5, h)), fragColor, fragCoord, 4.0);
    // }

    
    fragColor = vec4(c, 1.);

}