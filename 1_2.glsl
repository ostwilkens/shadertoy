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
    return min(1., (1. / distanceField) * 0.005) * 0.2;
}

float cube(vec3 p, float r) {
    return length(p - clamp(p, -r, r));
}

float sphere(vec3 p, float radius)
{
    return length(p) - radius;
}

float point(vec3 p, float r, float s) {
    return min(1., (1. / length(p)) * r) * s;
}

float smoothstep2(float a, float b, float v) {
    float op = (step(0., v) * 2.) - 1.;
    return smoothstep(a, b, abs(v)) * op;
}

vec4 scene(vec3 p)
{
    float n = iTime / BEAT;
    p.xy *= 1. + length(p.xy) * 0.2; // barrel distort
    p.xy /= 0.35; // zoom out
    // p.xy /= 0.2;

    n *= 0.1;

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
    p *= (0.95 + smoothstep(0., 0.2, mod(n - 0.43, 1.)) * 0.05); // pump
    float recoil = (0.9 + smoothstep(0., 0.9, mod(n - 0.2, 4.)) * 0.1);
    p *= recoil; // pump
    pR(p.xy, -n * 0.2 - recoil * 1.); // spin

    vec3 gp = p;

    // repeat
    vec3 cell = floor(p / 0.5);
    p.xy = mod(p.xy, 0.5);
    p.xy -= 0.25;
    float cellId = (1. - abs(cell.y * 3. - cell.x));
    pR(p.xy, -0.5 * cellId * PI);

    float orbInCell = abs((1. - abs(orbYPos - cell.y)) * (1. - abs(orbXPos - cell.x)));
    orbInCell *= cell.x == 0. ? 1. : 0. + cell.x == -1. ? 1. : 0.;
    orbInCell *= cell.y == 0. ? 1. : 0. + cell.y == -1. ? 1. : 0.;
    float orbNotInCell = 1. - orbInCell;
   
    p *= 1. + sin(((n + cellId * 4. + 7.) * PI * 0.5 + 0.95) / 4.) * 0.12; // slow zoom in<->out
    pR(p.xy, sin(n * 40.)*0.1 *stuckness * orbInCell); // spin
    pR(p.xy, (txnoise(p * 4. + n * 2.) - 0.5) * 1. * stuckness * orbInCell); // flare
    p *= (0.6 + smoothstep(0., 0.6, mod(wn - 0.58, 1.)) * 0.4) * orbInCell + orbNotInCell; // pump
    pR(p.yz, (p.x * 1. + n) * PI * 0.5); // x axis spiral
    
    vec3 p3 = p;

    float innerCube = cube(p, 0.048);
    float outerCube = cube(p, 0.05);
    float dCube = outerCube;
    float lCube = lightField(outerCube) - lightField(innerCube);

    vec3 p2 = gp;

    // point movement
    p2.x += warpnessX;
    p2.y += warpnessY;
    p2.x += sin(wn * 80.) * stuckness * 0.05;
    p2.y += sin(wn * 80.) * stuckness * 0.05;

    float point = point(p2, 0.006, 0.13);
    float dPoint = sphere(p2, 0.1);

    vec3 glow = vec3(0.0);
    glow += lCube * purple * 1.;
    glow += point * purple.grb;

    float glowStrength = (glow.r + glow.g + glow.b) * 0.7;
    float dist = 1. / (glowStrength * 10000.);

    return vec4(glow, dist);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (fragCoord.xy - (iResolution.xy * 0.5)) / iResolution.yy;
    uv *= 1. + length(uv) * 0.23; // barrel distort

    vec3 cameraOrigin = vec3(0., 0., -1.0);
    vec3 cameraTarget = cameraOrigin + vec3(0., 0., 1.);
    vec3 upDirection = vec3(0., 1.0, 0.);
    vec3 cameraDir = normalize(cameraTarget - cameraOrigin);

    vec3 cameraRight = normalize(cross(upDirection, cameraOrigin));
	vec3 cameraUp = cross(cameraDir, cameraRight);
    vec3 rayDir = normalize(cameraRight * uv.x + cameraUp * uv.y + cameraDir);


    vec3 c = vec3(0.0);
   
    vec3 pos = cameraOrigin;

    for(int i = 0; i < 200; i++)
    {
        vec4 ld = scene(pos);
        vec3 l = ld.rgb;
        float dist = ld.w;

        pos += rayDir * max(0.001, dist);
        c += l;

        if(dist > 10.)
            break;
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