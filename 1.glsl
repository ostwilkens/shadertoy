#iChannel0 "file://noise.png"
#iChannel0::WrapMode "Repeat"
#iChannel0::MagFilter "Linear"
#include "hg_sdf.glsl"

#define PRECISION 1000.
#define MAX_ITERATIONS 120
#define MAX_DIST 10.
#define EPSILON (1. / PRECISION + 0.0008)
#define time (iTime * (155. / 60.) / 1.)
const vec3 purple = normalize(vec3(0.298, 0.176, 0.459));

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

float shell(float d)
{
    return max(-d, d);
}

void scene2(vec3 p, float n, inout vec3 c, inout float d)
{
    n *= 0.5;
    p.z -= 0.5;
    pR(p.yz, 2.1);
    pR(p.xy, PI/2. + sin(n) * 1.5);
    add(d, c, max(-fDodecahedron(p, 0.15), fDodecahedron(p, 0.15)), purple, 5., 10.);
    pR(p.yx,  n);
    p.z += tan(n * PI) * 0.1;
    add(d, c, max(-fBox(p, vec3(0.2, 0.2, 0.2)), fBox(p, vec3(0.2, 0.2, 0.04))), purple.grb, 5., 0.4);
}

void scene3(vec3 p, float n, inout vec3 c, inout float d)
{
    n *= 0.5;

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
    
    add(d, c, shell(shapeC), purple, 6., 10.);
    add(d, c, shapeH, purple.grb, 9., 5.);
}

float s4exp(float off, float n)
{
    return smoothstep(off, 0., n) / n * 0.1 + 
        smoothstep(off, 0., mod(n, 6.)) * smoothstep(0., 6., mod(n, 6.)) * 0.1 + 
        0.2 / (10. + mod(n, 1.) * 40.);
}

void scene4(vec3 p, float n, inout vec3 c, inout float d)
{
    pR(p.zy, -0.2);
    pR(p.zx, n * 0.5);

    // p.x += tan(n) * 0.02;
    // pR(p.xy, n);
    add(d, c, fBoxCheap(p + vec3(0., 0., 0.1 + s4exp(1., n)), vec3(0.1, 0.1, 0.)), purple, 4., 20.);
    // p.y += tan(n) * 0.02;
    // pR(p.yz, n);
    add(d, c, fBoxCheap(p - vec3(0., 0., 0.1 + s4exp(2., n)), vec3(0.1, 0.1, 0.)), purple, 4., 20.);
    // p.z += tan(n) * 0.02;
    // pR(p.zx, n);
    add(d, c, fBoxCheap(p + vec3(0., 0.1 + s4exp(3., n), 0.), vec3(0.1, 0., 0.1)), purple, 4., 20.);
    // p.x += tan(n) * 0.02;
    // pR(p.xz, n);
    add(d, c, fBoxCheap(p - vec3(0., 0.1 + s4exp(4., n), 0.), vec3(0.1, 0., 0.1)), purple, 4., 20.);
    // p.y += tan(n) * 0.02;
    // pR(p.yx, n);
    add(d, c, fBoxCheap(p + vec3(0.1 + s4exp(5., n), 0., 0.), vec3(0., 0.1, 0.1)), purple, 4., 20.);
    // p.z += tan(n) * 0.02;
    // pR(p.zy, n);
    add(d, c, fBoxCheap(p - vec3(0.1 + s4exp(5.5, n), 0., 0.), vec3(0., 0.1, 0.1)), purple, 4., 20.);
}

void scene5(vec3 p, float n, inout vec3 c, inout float d)
{
    float b = 1. - smoothstep(0., 0.3, mod(n, 1.));
    float bi = (1. - b) + floor(n);

    p.z += n * 0.5;

    add(d, c, fBox2Cheap(p.xy + vec2(0., -0.3), vec2(0.2, 0.)), purple.grb, 1., 10.);

    p.y += sin(p.z * PI) * 0.01;
    pR(p.xy, p.z + 1.5);

    add(d, c, fBox2Cheap(p.xy + vec2(0., 0.1), vec2(0., 0.003 + sin(p.z) * 0.004)), purple.grb, 10., 1.);

    p.y += noise(p.xz * 8. + n * 1.5) * 0.05;
    // p.y += txnoise(vec3(p.xz * 8. + n * 1.5, 0.)) * 0.05;
    p.x -= noise(p.xz * 2. + n * 0.5 + 10.) * 0.1;
    p.y += noise(p.xz * 0.5 - n * 0.25) * 0.15;
    add(d, c, fBox2Cheap(p.xy + vec2(0., 0.15), vec2(0.1, 0.)), purple, 5., 5.);
}

void scene6(vec3 p, float n, inout vec3 c, inout float d)
{
    pR(p.yz, -0.2);
    vec3 p1 = p;

    pR(p.yz, 0.1 + p.z * 0.002);
    p.z += n * 1.;
    pR(p.xy, sin(p.z * 0.1) * 0.4);

    float f4 = noise(p.xz * 0.1) - 0.5;
    p.x += f4 * 4.;

    vec3 color = purple.grb * (1. - length(p1) * 0.008);

    add(d, c, shell(fBox2Cheap(p.xy + vec2(0., -100.05), vec2(0.05, 100.))), color, 2., 1.);
}

void scene7(vec3 p, float n, inout vec3 c, inout float d) {
    p.z -= sin(n * 0.3) * 0.5 + 0.2;

    vec3 p1 = p;
    pR(p1.xy, n * 0.1);
    add(d, c, fBox2(p1.xy + vec2(-0.0001, -0.001), vec2(0., 0.)), purple.grb, 10., 10.);
    // add(d, c, fBox2Cheap(p.xy + vec2(0., -0.25), vec2(0.25, 0.)), purple.grb, 2., 10.);


    // pR(p.yz, n);
    // add(d, c, fCircle(p, 0.3) - 0.01, purple, 2., 10.);

    pR(p.xy, n * 0.2 + sin(n) * 0.1);
    float c1 = pModPolar(p.xy, 6.);

    p.x -= 0.3 + sin(c1) * 0.1;
    p.z += tan(n + c1) * 0.2;
    pR(p.zy, PI / 2.);
    pR(p.xz, n * 0.6 + sin(n + c1) * 0.1 + c1);
    // add(d, c, shell(fSphere(p, 0.08)), purple, 1., 2.);
    // pR(p.yz, sin(n * 0.5));
    // pR(p.yz, -0.5);
    pModPolar(p.xz, 3.);
    pR(p.xy, tan(n + c1 + 0.3));
    add(d, c, shell(fTriPrism(p, vec2(0.06, 0.007))), purple, 8., 5.);

}

void scene8(vec3 p, float n, inout vec3 c, inout float d) {
    n *= 1.;
    n = mod(n, 122.);
    // n += 25.;
    n = 55.;

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

    add(d, c, sun, purple.brg, 15., 0.3);
    add(d, c, ocean, purple.grb, 2., 5.);
}

float point(vec3 p, float r, float s) {
    return min(1., (1. / length(p)) * r) * s;
}



// void scene9(vec3 p, float n, inout vec3 c, inout float d) {
//     add(d, c, glass(point(p, 10.1, 1.)), purple, 40., 1.);
// }

vec4 scene(vec3 p)
{
    float n = time;
    vec3 c = vec3(0.);
    float d = 1. / 0.;

    scene8(p, n, c, d);
    
    return vec4(c, d);
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
        // else if(i == MAX_ITERATIONS - 1)
        //     c = vec3(1., 0., 0.);
    }

    return c;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (fragCoord.xy - (iResolution.xy * 0.5)) / iResolution.yy;
    uv *= 1.5;
    vec3 origin = origin();
    vec3 upDirection = vec3(0., 1.0, 0.);
    vec3 cameraDir = normalize(target() - origin);
    vec3 cameraRight = normalize(cross(upDirection, origin));
	vec3 cameraUp = cross(cameraDir, cameraRight);
    vec3 rayDir = normalize(cameraRight * uv.x + cameraUp * uv.y + cameraDir);

    vec3 c = march(origin, rayDir);
    pp(c, uv);
   
    fragColor = vec4(c, 1.);
}