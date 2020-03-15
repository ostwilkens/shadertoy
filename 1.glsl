#iChannel0 "file://noise.png"
#iChannel0::WrapMode "Repeat"

#define PI 3.14159265358979323846
#define ITS 400.
#define BEAT 0.4477611940298507
const vec3 purple = normalize(vec3(0.298, 0.176, 0.459));

float txnoise(in vec3 x)
{
    // The MIT License
    // Copyright © 2013 Inigo Quilez
    // Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);
    vec2 uv = (p.xy+vec2(37.0,17.0)*p.z) + f.xy;
    vec2 rg = textureLod( iChannel0, (uv+0.5)/256.0, 0.0).yx;
    return mix( rg.x, rg.y, f.z );
}

void pR(inout vec2 p, float a) {
    // thx to hg (http://mercury.sexy/hg_sdf)
    p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

void pR45(inout vec2 p) {
    // thx to hg (http://mercury.sexy/hg_sdf)
	p = (p + vec2(p.y, -p.x))*sqrt(0.5);
}

float smin( float a, float b, float k )
{
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

float displacement(vec3 p, float k)
{
	return sin(k*p.x)*sin(k*p.y)*sin(k*p.z);
}

float lightField(float distanceField) {
    return min(1., (1. / distanceField) * 0.01) * 0.2;
}

float cube(vec3 p, float r) {
    vec3 field = clamp(p, -r, r);
    float cube = (1. / length(p - field)) * 0.01;
    return min(1., cube) * 0.2;
}

float sphere(vec3 p, float r) {
    return (1. - step(r, length(p) * 0.01)) * 0.01;
}

// float sphere(vec3 pos, float radius)
// {
//     float dist = length(pos) + radius;
//     return lightField(dist);
// }

float point(vec3 p, float r, float s) {
    return min(1., (1. / length(p)) * r) * s;
}

float smoothstep2(float a, float b, float v) {
    float op = (step(0., v) * 2.) - 1.;
    return smoothstep(a, b, abs(v)) * op;
}

vec3 scene(vec3 p)
{
    float n = iTime / BEAT;
    p.xy *= 1. + length(p.xy) * 0.2; // barrel distort
    p /= 0.9; // zoom out

    // n = mod(n, 4.); // loop dur
    // n += 4.; // start offset
    n += 53.1;
    // n *= 0.2;


    float wn = n * 0.2;
    float warpnessX = -0.25 + smoothstep(0.0, 0.6, mod(wn, 4.)) * 0.5 - smoothstep(2.0, 2.6, mod(wn, 4.)) * 0.5;
    float warpnessY = -0.25 + smoothstep(1.0, 1.6, mod(wn, 4.)) * 0.5 - smoothstep(3.0, 3.6, mod(wn, 4.)) * 0.5;
    float stuckness = (abs(warpnessX) + abs(warpnessY));
    stuckness = smoothstep(0.42, 0.5, stuckness) + 0.05;


    // p *= 1. + sin((n) / 2.) * 0.2; // slow zoom in<->out
    // p *= 1. + tan((n * PI + 0.95) / 4.) * 0.2; // warp in<->out
    // pR(p.xy, p.z * 1. + n * 0.05); // spin
    // pR(p.zx, txnoise(p * 2. + n)); // flare
    pR(p.xy, (txnoise(p * 4. + n) - 0.5) * 0.05 * stuckness); // flare
    // pR(p.yz, (p.x * 18. + n * 1.) * 0.1); // x axis spiral


    vec3 gp = p;
    
    pR(p.zx, sin(n * 0.25) * 0.3); // rot
    pR(p.zy, cos(n * 0.25) * 0.5); // rot
    // pR(p.zy, tan(n * 0.25) * 0.35); // rot
    // pR(p.xy, n * 0.1); // rot





    vec3 cell = floor(p);
    p.xy += 0.25;
    p *= 2.;
    p = fract(p + 0.5) - 0.5;
    p /= 2.;
    pR(p.xy, -0.5 * (1. - abs(cell.y * 3. - cell.x)) * PI);


    // pR(p.zx, n * 0.1); // rot
    // pR(p.xy, n * 0.1); // rot
   
    // p /= max(0.8, (1.-fract(n))) + 0.2; // pump
    
    // vec3 cell = floor(p);
    // p *= 2.;
    // p = fract(p + 0.5) - 0.5;
    // p /= 2.;
   
    // p *= smoothstep(0., 0.5, n); // initial warp in
    // p *= 1. + sin((n * PI + 0.95) / 4.) * 0.2; // slow zoom in<->out
    // p *= 1. + tan((n * PI + 0.95) / 4.) * 0.2; // warp in<->out
    // pR(p.xy, p.z * 2. + n); // spin
    // pR(p.zx, txnoise(p * 2. + n)); // flare
    // pR(p.xy, txnoise(p * 4. + n) - 0.5); // flare
    // pR(p.yz, (p.x * 4. + n * 1.) * 1.5); // x axis spiral
    // pR(p.xy, n); // rot
    // p.x /= 5.; // long x
    
    vec3 p3 = p;
    float rotatePipe = step(p3.y, 0.05) + step(p3.x, 0.05);
    pR(p3.xy, rotatePipe * PI * 0.5);
    float pipe = sphere(vec3(p3.x, 0., p3.z), 0.0002);
    pipe -= sphere(vec3(p3.x, 0., p3.z), 0.00015);

    float innerCube = cube(p, 0.048);
    float outerCube = cube(p, 0.05);
    outerCube = max(0., outerCube - pipe);
    float cube = outerCube - innerCube;

    vec3 p2 = gp;
    // p2.x += smoothstep(0., 1., abs(sin(n)));
    // float warpnessX = smoothstep2(0., 1., mod(n * 0.5 * step(1., n), 2.) -  mod((n + 2.) * 0.5, 2.)) * 0.25;

    p2.x += warpnessX;
    p2.y += warpnessY;
    p2.x += sin(wn * 100.) * stuckness * 0.05;
    p2.y += sin(wn * 100.) * stuckness * 0.05;

    // float warpness = tan(n * 0.5) * 0.03;
    // p2.y += warpness;
    // p2.x -= cos(n * 10.) * 0.05 * smoothstep(0.45, 0.5, 1. - abs(stuckness));
    float point = point(p2, 0.006, 0.15);
    

    pipe = max(0., pipe - point - innerCube * 0.1);
    cube = max(0., cube - pipe);
    

    vec3 glow = vec3(0.0);
    glow += cube * purple * 1.8;
    glow += pipe * purple * 6.;
    glow += point * purple.grb;

    vec3 c = glow;
    c *= 200. / ITS;
    // c *= smoothstep(0., 3., n); // fade in
    return c;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (fragCoord.xy - (iResolution.xy * 0.5)) / iResolution.yy;
    uv *= 1. + length(uv) * 0.2; // barrel distort

    vec3 c = vec3(0.0);
   
    float hr = 0.2;
    for(float z = -hr; z < hr; z += 1.0 / ITS)
    {
        c += scene(vec3(uv, z));
    }
    c = clamp(c, 0., 1.);
    
    // c += sin(iTime * 0.4) * 0.08; // ambience

	// grading
    c = clamp(c, 0.02, 1.);
	c -= 0.02;
	c *= 1.1;
	c = sqrt(c);
	c = c*c*(2.5-1.5*c*c); // contrast
	c = pow(c, vec3(1.0,0.96,1.0)); // soft green
	c *= vec3(1.08,0.99,0.99); // tint red
	c.z = (c.z+0.05)/1.05; // bias blue
	c = mix(c, c.yyy, 0.12); // desaturate
	c = 0.06 + (c * 0.9);
    
    c -= smoothstep(0.55, 1.3, abs(uv.x)) * 0.2; // vignette x
    c -= smoothstep(0.17, 0.7, abs(uv.y)) * 0.2; // vignette y
    c -= step(0.515, abs(uv.x)); // letterbox x
    c -= step(0.35, abs(uv.y)); // letterbox y
    c = clamp(c, 0., 1.);
    c = max(c, txnoise(vec3(uv.xy, 0.) * 700.) * 0.02 + 0.005); // texturize black
    c += txnoise(vec3(uv.xy, sin(cos(iTime) * 1000.)) * 1000.) * 0.04; // noise
   
    fragColor = vec4(c, 1.);
}