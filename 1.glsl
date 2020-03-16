#iChannel0 "file://noise.png"
#iChannel0::WrapMode "Repeat"

#define PI 3.14159265358979323846
#define ITS 150
#define BEAT 0.4477611940298507
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

void pR(inout vec2 p, float a) {
    // thx to hg (http://mercury.sexy/hg_sdf)
    p = cos(a)*p + sin(a)*vec2(p.y, -p.x);
}

void pR45(inout vec2 p) {
    // thx to hg (http://mercury.sexy/hg_sdf)
	p = (p + vec2(p.y, -p.x))*sqrt(0.5);
}

float smin( float a, float b, float dist )
{
    float h = clamp( 0.5+0.5*(b-a)/dist, 0.0, 1.0 );
    return mix( b, a, h ) - dist*h*(1.0-h);
}

float displacement(vec3 p, float dist)
{
	return sin(dist*p.x)*sin(dist*p.y)*sin(dist*p.z);
}

float lightField(float distanceField) {
    return min(1., (1. / distanceField) * 0.01) * 0.2;
}

float cube(vec3 p, float r) {
    vec3 field = clamp(p, -r, r);
    float cube = (1. / length(p - field)) * 0.0017;
    return min(1., cube) * 0.2;
}

float sphere(vec3 p, float r) {
    return (1. - step(r, length(p) * 0.01)) * 0.01;
}

float point(vec3 p, float r, float s) {
    return min(1., (1. / length(p)) * r) * s;
}

float smoothstep2(float a, float b, float v) {
    float op = (step(0., v) * 2.) - 1.;
    return smoothstep(a, b, abs(v)) * op;
}


void worldR(inout vec3 p, float op) {
    float n = iTime / BEAT;
    // pR(p.yx, 1.0 * op);
    pR(p.zx, sin(n * 0.25) * 0.3 * op); // rot
    pR(p.zy, cos(n * 0.25) * 0.2 * op); // rot
}

vec3 scene(vec3 p)
{
    float n = iTime / BEAT;
    worldR(p, 1.);
    p.xy *= 1. + length(p.xy) * 0.2; // barrel distort
    p /= 0.35; // zoom out
    // p.xy /= 0.3;


    // n = mod(n, 4.); // loop dur
    // n += 4.; // start offset
    // n += 53.1;
    // n *= 0.2;

    // point movement calc
    float wn = n / 4. + 0.46;
    float warpnessX = -0.25 + smoothstep(0.45, 0.6, mod(wn, 4.)) * 0.5 - smoothstep(2.45, 2.6, mod(wn, 4.)) * 0.5;
    float warpnessY = -0.25 + smoothstep(1.45, 1.6, mod(wn, 4.)) * 0.5 - smoothstep(3.45, 3.6, mod(wn, 4.)) * 0.5;
    float orbXPos = step(0., warpnessX);
    float orbYPos = step(0., warpnessY);
    float stuckness = (abs(warpnessX) + abs(warpnessY));
    stuckness = smoothstep(0.42, 0.5, stuckness) + 0.05;


    // global modifiers
    p *= 1. + sin(n*0.15) * 0.3; // slow zoom in<->out
    // p *= 1. + tan((n * PI * 0.5) / 8.) * 0.02; // warp in<->out
    // pR(p.zx, n * 0.1); // spin
    // pR(p.zx, txnoise(p * 2. + n)); // flare
    // pR(p.xy, (txnoise(p * 4. + n) - 0.5) * 0.15 * stuckness); // flare
    // pR(p.yz, (p.x * 18. + n * 1.) * 0.1); // x axis spiral
    p *= (0.95 + smoothstep(0., 0.2, mod(n - 0.43, 1.)) * 0.05); // pump
    float recoil = (0.9 + smoothstep(0., 0.7, mod(n - 0.2, 4.)) * 0.1);
    p *= recoil; // pump
    pR(p.xy, -n * 0.2 + recoil * 2.); // spin

    // p.xy -= vec2(warpnessX, warpnessY) * 1.; // follow cam
    

    // p.xy += 2.;

    // p.xy += 0.4;
    // // repeat orb
    // p.xy = mod(p.xy, 1.);
    // p.xy -= 0.5;

    vec3 gp = p;

    
    // pR(p.zx, sin(n * 0.25) * 0.3); // rot
    // pR(p.zy, cos(n * 0.25) * 0.5); // rot
    // // pR(p.zy, tan(n * 0.25) * 0.35); // rot
    // // pR(p.xy, n * 0.1); // rot



    // repeat
    vec3 cell = floor(p / 0.5);
    p.xy = mod(p.xy, 0.5);
    p.xy -= 0.25;
    float cellId = (1. - abs(cell.y * 3. - cell.x));
    pR(p.xy, -0.5 * cellId * PI);


    float orbInCell = abs((1. - abs(orbYPos - cell.y)) * (1. - abs(orbXPos - cell.x)));
    // orbInCell = clamp(orbInCell, 0., 1.);
    orbInCell *= cell.x == 0. ? 1. : 0. + cell.x == -1. ? 1. : 0.;
    orbInCell *= cell.y == 0. ? 1. : 0. + cell.y == -1. ? 1. : 0.;
    // orbInCell = clamp(orbInCell, 0., 1.);
    float orbNotInCell = 1. - orbInCell;
    // pR(p.xy, orbInCell);
    // p.z += orbInCell;


    // pR(p.zx, (n + abs(cell.y * 3. - cell.x) + 3.) * 0.5); // rot
    // pR(p.xy, 0.5 * PI * n); // rot
   
    // p /= max(0.8, (1.-fract(n))) + 0.2; // pump
    
    // vec3 cell = floor(p);
    // p *= 2.;
    // p = fract(p + 0.5) - 0.5;
    // p /= 2.;
   
    // p *= smoothstep(0., 0.5, n); // initial warp in
    p *= 1. + sin(((n + cellId * 4. + 7.) * PI * 0.5 + 0.95) / 4.) * 0.12; // slow zoom in<->out
    // p *= 1. + tan((n * PI + 0.95) / 4.) * 0.2; // warp in<->out
    // pR(p.xy, p.z * 2. + 0.5 * PI * (n + 3.)); // spin
    pR(p.xy, sin(n * 40.)*0.1 *stuckness * orbInCell); // spin
    // pR(p.zx, txnoise(p * 2. + n)); // flare
    // pR(p.xy, txnoise(p * 4. + n) - 0.5); // flare
    pR(p.xy, (txnoise(p * 4. + n * 2.) - 0.5) * 1. * stuckness * orbInCell); // flare
    p *= (0.6 + smoothstep(0., 0.6, mod(wn - 0.58, 1.)) * 0.4) * orbInCell + orbNotInCell; // pump
    // p *= (0.7 + smoothstep(0., 0.4, mod(n - 0.1, 3.9)) * 0.35) * orbInCell + orbNotInCell;
    pR(p.yz, (p.x * 1. + n) * PI * 0.5); // x axis spiral
    // pR(p.xy, 0.5 * PI * n); // rot
    
    vec3 p3 = p;
    float rotatePipe = step(p3.y, 0.05) + step(p3.x, 0.05);
    pR(p3.xy, rotatePipe * PI * 0.5);
    float pipe = sphere(vec3(p3.x, 0., p3.z), 0.00022);
    pipe -= sphere(vec3(p3.x, 0., p3.z), 0.00014);

    float innerCube = cube(p, 0.048);
    float outerCube = cube(p, 0.05);
    float cube = outerCube - innerCube;

    vec3 p2 = gp;

    // point movement
    p2.x += warpnessX;
    p2.y += warpnessY;
    p2.x += sin(wn * 80.) * stuckness * 0.05;
    p2.y += sin(wn * 80.) * stuckness * 0.05;

    float point = point(p2, 0.006, 0.13);
    
    pipe = max(0., pipe - innerCube * 0.1 - step(0.12, length(p)));
    cube = max(0., cube - pipe);

    vec3 glow = vec3(0.0);
    glow += cube * purple * 1.4;
    glow += pipe * purple * 6.;
    glow += point * purple.grb;

    vec3 c = glow;
    c *= 0.55;
    // c *= 200. / float(ITS);
    // c *= smoothstep(0., 3., n); // fade in
    return c;
}


void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (fragCoord.xy - (iResolution.xy * 0.5)) / iResolution.yy;
    uv *= 1. + length(uv) * 0.23; // barrel distort

    vec3 cameraOrigin = vec3(0., 0., -0.1);
    vec3 cameraTarget = cameraOrigin + vec3(0., 0., 1.);
    vec3 upDirection = vec3(0., 1.0, 0.);
    vec3 cameraDir = normalize(cameraTarget - cameraOrigin);

    worldR(cameraDir, -1.);

    vec3 cameraRight = normalize(cross(upDirection, cameraOrigin));
	vec3 cameraUp = cross(cameraDir, cameraRight);
    vec3 rayDir = normalize(cameraRight * uv.x + cameraUp * uv.y + cameraDir);


    vec3 c = vec3(0.0);
   
    vec3 pos = cameraOrigin;

    float dist = 0.6;

    pos.z -= dist;
    worldR(pos, -1.);
    
    pos += rayDir * (dist + 0.075);

    for(int i = 0; i < ITS; i++)
    {
        pos += rayDir * (0.15 / float(ITS));

        c += scene(pos);
    }
    c = clamp(c, 0., 1.);
    
	// grading
    c = clamp(c, 0.02, 1.);
	c -= 0.02;
	c *= 1.2;
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
    c = max(c, txnoise(vec3(uv.xy, 0.) * 700.) * 0.025 + 0.01); // texturize black
    c += txnoise(vec3(uv.xy, sin(cos(iTime) * 1000.)) * 1000.) * 0.04; // noise
   
    fragColor = vec4(c, 1.);
}