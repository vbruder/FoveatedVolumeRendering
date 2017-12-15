#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

#define ERT_THRESHOLD 0.98

constant sampler_t linearSmp = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_CLAMP |
                                CLK_FILTER_LINEAR;
constant sampler_t nearestSmp = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_CLAMP |
                                CLK_FILTER_NEAREST;

// intersect ray with a box
// http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm
int intersectBox(float3 rayOrig, float3 rayDir, float *tnear, float *tfar)
{
    // compute intersection of ray with all six bbox planes
    float3 invRay = native_divide((float3)(1.0f), rayDir);
    float3 tBot = invRay * ((float3)(-1.0f) - rayOrig);
    float3 tTop = invRay * ((float3)(1.0f) - rayOrig);

    // re-order intersections to find smallest and largest on each axis
    float3 tMin = min(tTop, tBot);
    float3 tMax = max(tTop, tBot);

    // find the largest tMin and the smallest tMax
    float maxTmin = max(max(tMin.x, tMin.y), max(tMin.x, tMin.z));
    float minTmax = min(min(tMax.x, tMax.y), min(tMax.x, tMax.z));

    *tnear = maxTmin;
    *tfar = minTmax;

    return (int)(minTmax > maxTmin);
}

int intersectBBox(float3 rayOrig, float3 rayDir, float3 lower, float3 upper, float *tnear, float *tfar)
{
    // compute intersection of ray with all six bbox planes
    float3 invRay = native_divide((float3)(1.0f), rayDir);
    float3 tBot = invRay * (lower - rayOrig);
    float3 tTop = invRay * (upper - rayOrig);

    // re-order intersections to find smallest and largest on each axis
    float3 tMin = min(tTop, tBot);
    float3 tMax = max(tTop, tBot);

    // find the largest tMin and the smallest tMax
    float maxTmin = max(max(tMin.x, tMin.y), max(tMin.x, tMin.z));
    float minTmax = min(min(tMax.x, tMax.y), min(tMax.x, tMax.z));

    *tnear = maxTmin;
    *tfar = minTmax;

    return (int)(minTmax > maxTmin);
}


int intersectPlane(const float3 rayOrigin, const float3 rayDir,
                   const float3 planeNormal, const float3 planePos, float *t)
{
    float denom = dot(rayDir, planeNormal);
    if (denom > 1e-6)   // check if parallel
    {
        float3 origin2point = planePos - rayOrigin;
        *t = dot(origin2point, planeNormal) / denom;
        return (t >= 0);
    }
    return false;
}

// Compute gradient using central difference: f' = ( f(x+h)-f(x-h) ) / 2*h
float3 gradientCentralDiff(read_only image3d_t vol, const float4 pos)
{
    float3 volResf = convert_float3(get_image_dim(vol).xyz);
    float3 offset = native_divide((float3)(1.0f, 1.0f, 1.0f), volResf);
    float3 s1;
    float3 s2;
    s1.x = read_imagef(vol, nearestSmp, pos + (float4)(-offset.x, 0, 0, 0)).x;
    s2.x = read_imagef(vol, nearestSmp, pos + (float4)(+offset.x, 0, 0, 0)).x;
    s1.y = read_imagef(vol, nearestSmp, pos + (float4)(0, -offset.y, 0, 0)).x;
    s2.y = read_imagef(vol, nearestSmp, pos + (float4)(0, +offset.y, 0, 0)).x;
    s1.z = read_imagef(vol, nearestSmp, pos + (float4)(0, 0, -offset.z, 0)).x;
    s2.z = read_imagef(vol, nearestSmp, pos + (float4)(0, 0, +offset.z, 0)).x;
    return (s2 - s1) / (2.f * offset);
}

// specular part of blinn-phong shading model
float3 specularBlinnPhong(float3 lightColor, float specularExp, float3 materialColor,
                          float3 normal, float3 toLightDir, float3 toCameraDir)
{
    float3 h = toCameraDir + toLightDir;

    // check for special case where the light source is exactly opposite
    // to the view direction, i.e. the length of the halfway vector is zero
    if (dot(h, h) < 1.e-6f) // check for squared length
        return (float3)(0.0f);

    h = normalize(h);
    return materialColor * lightColor * native_powr(max(dot(normal, h), 0.f), specularExp);
}

// simple illumination based on central differences
float3 illumination(read_only image3d_t vol, const float4 pos, float3 diffuse,float3 toCameraDir)
{
    float3 n = fast_normalize(gradientCentralDiff(vol, pos));
    float3 l = fast_normalize((float3)(20.0f, 100.0f, 20.0f) - pos.xyz);

    float3 amb = diffuse;
    float3 diff = diffuse * max(0.f, dot(n, l));
    float3 spec = specularBlinnPhong((float3)(1.f), 100.f, (float3)(1.f), n, l, toCameraDir);

    return (amb + diff + spec) * 0.5f;
}

// check for border edges
bool checkBoundingBox(float3 pos, float3 voxLen, float2 bound)
{
    if (
        (pos.x < voxLen.x     && pos.z < bound.x + voxLen.z) ||
        (pos.x < voxLen.x     && pos.y < voxLen.y) ||
        (pos.y < voxLen.y     && pos.z < bound.x + voxLen.z) ||
        (pos.x > 1.f-voxLen.x && pos.z < bound.x + voxLen.z) ||
        (pos.y > 1.f-voxLen.y && pos.z < bound.x + voxLen.z) ||
        (pos.x > 1.f-voxLen.x && pos.z > bound.y - voxLen.z) ||
        (pos.y > 1.f-voxLen.y && pos.z > bound.y - voxLen.z) ||
        (pos.x < voxLen.x     && pos.z > bound.y - voxLen.z) ||
        (pos.y < voxLen.y     && pos.z > bound.y - voxLen.z) ||
        (pos.x > 1.f-voxLen.x && pos.y < voxLen.y) ||
        (pos.x > 1.f-voxLen.x && pos.y > 1.f -voxLen.y) ||
        (pos.x < voxLen.x     && pos.y > 1.f -voxLen.y)
       )
        return true;
    else
        return false;
}

/**
 * direct volume raycasting kernel
 */
__kernel void volumeRender(  __read_only image3d_t volData
                           , __read_only image1d_t tffData     // constant transfer function values
                           , __read_only image3d_t volBrickData
                           , __write_only image2d_t outData
                           , const float samplingRate
                           , const float16 viewMat
                           , const uint orthoCam
                           , const uint useIllum
                           , const uint useBox
                           , const uint useLinear
                           , const float4 background
                           , __read_only image1d_t tffPrefix
                           )
{
    int2 globalId = (int2)(get_global_id(0), get_global_id(1));
    if(any(globalId >= get_image_dim(outData)))
        return;

    // pseudo random number [0,1] for ray offsets to avoid moire patterns
    float iptr;
    float rand = 0.1f*fract(sin(dot(convert_float2(globalId),
                       (float2)(12.9898f, 78.233f))) * 43758.5453f, &iptr);

    float maxRes = (float)max(get_image_dim(volData).x, get_image_dim(volData).z);
//    float stepSize = native_divide(2.f, maxRes*samplingRate); // normalization

    int2 texCoords = globalId;
    float aspectRatio = native_divide((float)get_global_size(1), (float)(get_global_size(0)));
    aspectRatio = min(aspectRatio, native_divide((float)get_global_size(0), (float)(get_global_size(1))));
    int maxSize = max(get_global_size(0), get_global_size(1));
    float2 imgCoords;
    imgCoords.x = native_divide((globalId.x + 0.5f), convert_float(maxSize)) * 2.f;
    imgCoords.y = native_divide((globalId.y + 0.5f), convert_float(maxSize)) * 2.f;
    // calculate correct offset based on aspect ratio
    imgCoords -= get_global_size(0) > get_global_size(1) ?
                        (float2)(1.0f, aspectRatio) : (float2)(aspectRatio, 1.0);
    imgCoords.y *= -1.f;   // flip y coord

    float4 rayDir = (float4)(0.f);
    
    // z position of view plane is -1.0 to fit the cube to the screen quad when axes are aligned, 
    // zoom is -1 and the data set is uniform in each dimension
    // (with FoV of 90° and near plane in range [-1,+1]).
    float4 nearPlanePos = fast_normalize((float4)(imgCoords, -1.0f, 0.0f));
    // transform nearPlane from view space to world space
    rayDir.x = dot(viewMat.s0123, nearPlanePos);
    rayDir.y = dot(viewMat.s4567, nearPlanePos);
    rayDir.z = dot(viewMat.s89ab, nearPlanePos);
    rayDir.w = dot(viewMat.scdef, nearPlanePos);
    
    // camera position in world space (ray origin) is translation vector of view matrix
    float4 camPos = (float4)(viewMat.s37b, 1.0f);
    rayDir = fast_normalize(rayDir);

    if (orthoCam)
    {
        camPos = (float4)(viewMat.s37b, 1.0f);
        float4 viewPlane_x = viewMat.s048c;
        float4 viewPlane_y = viewMat.s159d;
        float4 viewPlane_z = viewMat.s26ae;
        rayDir = -viewPlane_z;
        rayDir.w = 0;
        rayDir = fast_normalize(rayDir);
        nearPlanePos = camPos + imgCoords.x*viewPlane_x + imgCoords.y*viewPlane_y;
        nearPlanePos *= length(camPos);
        camPos = nearPlanePos;
    }

    float tnear = FLT_MIN;
    float tfar = FLT_MAX;
    int hit = 0;
    // bbox from (-1,-1,-1) to (+1,+1,+1)
    hit = intersectBox(camPos.xyz, rayDir.xyz, &tnear, &tfar);
    if (!hit || tfar < 0)
    {
        write_imagef(outData, texCoords, background);
        return;
    }

    float sampleDist = tfar - tnear;
    if (sampleDist <= 0.f)
        return;
    float stepSize = min(sampleDist, sampleDist / (samplingRate*length(sampleDist*rayDir.xyz*convert_float3(get_image_dim(volData).xyz))));
    float samples = ceil(sampleDist/stepSize);
    stepSize = sampleDist/samples;

    // raycast parameters
    tnear = max(0.f, tnear);    // clamp to near plane
    float4 result = background;
    float alpha = 0.f;
    float4 pos = (float4)(0);
    float density = 0.f;
    float4 tfColor = (float4)(0);
    float opacity = 0.f;
    float t = tnear + rand*stepSize;    // offset by 'random' distance to avoid moiré pattern

    int3 volRes = get_image_dim(volData).xyz;
    float3 voxLen = (float3)(1.f) / convert_float3(volRes);
    int3 bricksRes = get_image_dim(volBrickData).xyz;
    float3 brickLen = (float3)(1.f) / convert_float3(bricksRes);
    float refSamplingInterval = 1.f / samplingRate;
    float precisionDiv = 1.f;
    if (get_image_channel_data_type(volData) == CLK_UNORM_INT16)
        precisionDiv = 8.f;

    // DDA initialization
    float3 invRay = 1.f/rayDir.xyz;
    int3 step = convert_int3(sign(rayDir.xyz));
    if (rayDir.x == 0.f)
    {
        invRay.x = FLT_MAX;
        step.x = 1;
    }
    if (rayDir.y == 0.f)
    {
        invRay.y = FLT_MAX;
        step.y = 1;
    }
    if (rayDir.z == 0.f)
    {
        invRay.z = FLT_MAX;
        step.z = 1;
    }

    float3 deltaT = convert_float3(step)*(brickLen*2.f*invRay);
    float3 voxIncr = (float3)0;

    // convert ray starting point to cell coordinates
    float3 rayOrigCell = (camPos.xyz + rayDir.xyz * tnear) - (float3)(-1.f);
    int3 cell = clamp(convert_int3(floor(rayOrigCell / (2.f*brickLen))),
                        (int3)(0), convert_int3(bricksRes.xyz) - 1);

    // add +1 to cells if ray dir component is negative: rayDir >= 0 ? (-1) : 0
    float3 tv = tnear + (convert_float(cell - isgreaterequal(rayDir.xyz, (float3)(0))) * (2.f*brickLen) - rayOrigCell) * invRay;
    int3 exit = step * bricksRes.xyz;
    if (exit.x < 0) exit.x = -1;
    if (exit.y < 0) exit.y = -1;
    if (exit.z < 0) exit.z = -1;

    uint i = 0;
    float t_exit = 0;
    float t_entry = 0;
    int cnt = 0;
    // raycasting loop: front to back raycasting with early ray termination
    while (t <= tfar)
    {
        float2 minMaxDensity = read_imagef(volBrickData, nearestSmp, (int4)(cell, 0)).xy;

        // increment to next brick
        voxIncr.x = (tv.x <= tv.y) && (tv.x <= tv.z) ? 1 : 0;
        voxIncr.y = (tv.y <= tv.x) && (tv.y <= tv.z) ? 1 : 0;
        voxIncr.z = (tv.z <= tv.x) && (tv.z <= tv.y) ? 1 : 0;
        cell += convert_int3(voxIncr) * step;    // [0; res-1]

        t_exit = (1.f*dot((float3)(1), tv * voxIncr));
        if (t_exit < t)
        {
//            write_imagef(outData, texCoords, (float4)(1,0,0,1));
//            return;
            //t_exit = t + 0.9*stepSize;
            t = t_exit - 1.1f*stepSize;
            ++cnt;
        }
        tv += voxIncr*deltaT;

        // skip bricks that contain only fully transparent voxels
        float alphaMax = read_imagef(tffData, nearestSmp, minMaxDensity.y).w;
        if (alphaMax < 0.001f)
        {
            uint prefixMin = read_imageui(tffPrefix, nearestSmp, minMaxDensity.x).x;
            uint prefixMax = read_imageui(tffPrefix, nearestSmp, minMaxDensity.y).x;
            if (prefixMin == prefixMax)
            {
                t = t_exit;
                continue;
            }
        }
        while (t < t_exit)
        {
            // t += tnear + stepSize*i;     // recalculate t to avoid numerical drift
            pos = camPos + t*rayDir;
            pos = pos * 0.5f + 0.5f;    // normalize to [0,1]

            density = useLinear ? read_imagef(volData, linearSmp, pos).x :
                                  read_imagef(volData, nearestSmp, pos).x;
            density = clamp(density, 0.f, 1.f);
            tfColor = read_imagef(tffData, linearSmp, density);  // map density to color
            if (useIllum)
                tfColor.xyz = illumination(volData, pos, tfColor.xyz, fast_normalize(camPos.xyz - pos.xyz));
            tfColor.xyz = background.xyz - tfColor.xyz;

            // Taylor expansion approximation
            opacity = 1.f - native_powr(1.f - tfColor.w, refSamplingInterval);
            result.xyz = result.xyz - tfColor.xyz * opacity * (1.f - alpha);
            alpha = alpha + opacity * (1.f - alpha);

            if (alpha > ERT_THRESHOLD || t >= tfar) break;   // early ray termination check
            ++i;
            t += stepSize;
        }

        if (alpha > ERT_THRESHOLD) break;
        if (any(cell == exit)) break;
        t = t_exit;
    }

    // draw bounding box
    if (useBox && checkBoundingBox(pos.xyz, voxLen, (float2)(0.f, 1.f)))
    {
        result.xyz = fabs((float3)(1.f) - background.xyz);
        opacity = 1.f;
    }

    result.w = alpha;
    write_imagef(outData, texCoords, result);
}

#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

//************************** Generate brick volume ***************************

__kernel void generateBricks(  __read_only image3d_t volData
                             , __read_only image1d_t tffData     // constant transfer function values
                             , __write_only image3d_t volBrickData
                            )
{
    int3 coord = (int3)(get_global_id(0), get_global_id(1), get_global_id(2));
    if(any(coord >= get_image_dim(volBrickData).xyz))
        return;

    int3 voxPerCell = convert_int3(ceil(convert_float4(get_image_dim(volData))/
                                          convert_float4(get_image_dim(volBrickData))).xyz);
    int3 volCoordLower = voxPerCell * coord;
    int3 volCoordUpper = clamp(volCoordLower + voxPerCell, (int3)(0), get_image_dim(volData).xyz);

    float maxVal = 0.f;
    float minVal = FLT_MAX;
    float value = 0.f;

    for (int k = volCoordLower.z; k < volCoordUpper.z; ++k)
    {
        for (int j = volCoordLower.y; j < volCoordUpper.y; ++j)
        {
            for (int i = volCoordLower.x; i < volCoordUpper.x; ++i)
            {
                value = read_imagef(volData, nearestSmp, (int4)(i, j, k, 0)).x;  // [0; 1]
                minVal = min(minVal, value);
                maxVal = max(maxVal, value);
            }
        }
    }

    write_imagef(volBrickData, (int4)(coord, 0), (float4)(minVal, maxVal, 0, 0));
}
