#iChannel0 "file://noise.png"
#iChannel0::WrapMode "Repeat"
#include "hg_sdf.glsl"

#define SHARPNESS 100.
#define PRECISION 100.
#define MAX_ITERATIONS 100
#define MAX_DIST 10.
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

float glass(float x) {
    return (1. / max(0., x)) / SHARPNESS;
}

float lightness(float df, float density)
{
    return min(1., glass(df)) * density * 10.;
}

void render(inout float d, inout vec3 c, float d2, vec3 color, float density)
{
    d = min(d, d2);
    c += lightness(d2, density) * color;
}

void render(inout float d, inout vec3 c, float d2, vec3 color)
{
    render(d, c, d2, color, 1.);
}

vec4 scene(vec3 p)
{
    float n = iTime * 0.1 + 4.;
    n = 4.;
    float d = 1. / 0.;
    vec3 c = vec3(0.);

    p.z -= 0.5;

    pR(p.xy, + n);
    pR(p.yz, + n);
    pR(p.zx, + n);

    render(d, c, fBox(p, vec3(0.1)), purple);
    // p.xy += 0.2;
    // render(d, c, fBox(p, vec3(0.1)), purple);
    // p.xy -= 0.4;
    // render(d, c, fBox(p, vec3(0.1)), purple.gbr);

    return vec4(c, d);
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

    for(int i = 0; i < MAX_ITERATIONS; i++)
    {
        vec4 s = scene(pos);
        float d = max(s.w, (1. + (0.1 * pos.z)) / PRECISION);

        pos += rayDir * d;
        c += s.rgb / (20. + PRECISION);

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
    c = clamp(c, 0., 1.);
    c = max(c, txnoise(vec3(uv.xy, 0.) * 700.) * 0.02); // texturize black
    c += txnoise(vec3(uv.xy, sin(cos(iTime) * 1000.)) * 1000.) * 0.04; // noise

    return c;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    vec2 uv = (fragCoord.xy - (iResolution.xy * 0.5)) / iResolution.yy;

    vec3 c = march(uv);
    c = post(c, uv);
   
    fragColor = vec4(c, 1.);
}