#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

#define ERT_THRESHOLD 0.98

constant sampler_t linearSmp = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_CLAMP |
                                CLK_FILTER_LINEAR;
constant sampler_t nearestSmp = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_CLAMP |
                                CLK_FILTER_NEAREST;

// Lambert shading
float4 shading(const float3 n, const float3 l, const float3 v)
{
    float4 ambient = (float4)(1.0f, 1.0f, 1.0f, 0.2f);
    float4 diffuse = (float4)(1.0f, 1.0f, 1.0f, 0.8f);

    float intensity = fabs(dot(n, l)) * diffuse.w;
    diffuse *= intensity;

    float4 color4 = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    color4.xyz = ambient.xyz * ambient.w + diffuse.xyz;
    return color4;
}


// intersect ray with a box
// http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm
int intersectBox(float4 rayOrig, float4 rayDir, float *tnear, float *tfar)
{
    // compute intersection of ray with all six bbox planes
    float4 invRay = native_divide((float4)(1.0f, 1.0f, 1.0f, 1.0f), rayDir);
    float4 tBot = invRay * ((float4)(-1.0f, -1.0f, -1.0f, 1.0f) - rayOrig);
    float4 tTop = invRay * ((float4)(1.0f, 1.0f, 1.0f, 1.0f) - rayOrig);

    // re-order intersections to find smallest and largest on each axis
    float4 tMin = min(tTop, tBot);
    float4 tMax = max(tTop, tBot);

    // find the largest tMin and the smallest tMax
    float maxTmin = max(max(tMin.x, tMin.y), max(tMin.x, tMin.z));
    float minTmax = min(min(tMax.x, tMax.y), min(tMax.x, tMax.z));

    *tnear = maxTmin;
    *tfar = minTmax;

    return (int)(minTmax > maxTmin);
}


float illumination(image3d_t vol, const float4 pos)
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
    float3 n = fast_normalize(s2 - s1).xyz;
    float3 l = fast_normalize((float4)(2.0f, 10.0f, 2.0f, 0.0f) - pos).xyz;

    return fabs(dot(n, l));
    //return shading(n, (normalize((float4)(1.0f, -1.0f, -1.0f, 0.0f) - pos)).xyz, (pos - rayDir));
}


uint getui4(uint4 v, int id)
{
    if (id == 0) return v.x;
    if (id == 1) return v.y;
    if (id == 2) return v.z;
    if (id == 3) return v.w;
}

float getf4(float4 v, int id)
{
    if (id == 0) return v.x;
    if (id == 1) return v.y;
    if (id == 2) return v.z;
    if (id == 3) return v.w;
}

float getf8(float8 v, int id)
{
    if (id == 0) return v.s0;
    if (id == 1) return v.s1;
    if (id == 2) return v.s2;
    if (id == 3) return v.s3;
    if (id == 4) return v.s4;
    if (id == 5) return v.s5;
    if (id == 6) return v.s6;
    if (id == 7) return v.s7;
}

bool approxEq(float x, float y)
{
    float delta = 0.005f;
    return (x < y + delta) && (y < x + delta);
}

bool approxEq2(float2 a, float2 b)
{
    return approxEq(a.x, b.x) && approxEq(a.y, b.y);
}


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


bool checkEdges(float3 pos, float3 boxMin, float3 boxMax)
{
    return
       approxEq2(pos.xz, boxMin.xz)
    || approxEq2(pos.yz, boxMin.yz)
    || approxEq2(pos.xz, boxMax.xz)
    || approxEq2(pos.yz, boxMax.yz)

    || approxEq2(pos.xz, (float2)(boxMin.x, boxMax.z))
    || approxEq2(pos.yz, (float2)(boxMin.y, boxMax.z))
    || approxEq2(pos.xz, (float2)(boxMax.x, boxMin.z))
    || approxEq2(pos.yz, (float2)(boxMax.y, boxMin.z));
}


/**
 * direct volume raycasting kernel
 */
__kernel void volumeRender(__read_only image3d_t volData,
                           __write_only image2d_t outData,
                           __read_only image1d_t tffData,     // constant transfer function values
                           const float stepSizeFactor,
                           __constant float *viewMat,
                           const uint orthoCam,
                           const uint useIllum,
                           const uint useBox,
                           const uint useLinear,
                           const float4 background
                           )
{
    int2 globalId = (int2)(get_global_id(0), get_global_id(1));
    if(any(globalId >= get_image_dim(outData)))
        return;

    float maxRes = (float)max(get_image_dim(volData).x, get_image_dim(volData).z);
    float stepSize = native_divide(stepSizeFactor, maxRes);
    stepSize *= 8.0f; // normalization to octile

    int2 texCoords = globalId;
    float2 imgCoords;
    imgCoords.x = native_divide((globalId.x + 0.5f), (float)(get_global_size(0))) * 2.0f - 1.0f;
    imgCoords.y = native_divide((globalId.y + 0.5f), (float)(get_global_size(1))) * 2.0f - 1.0f;
    imgCoords.y *= -1.0f;

    float4 rayDir = (float4)(0);
    float tnear = 0.0f;
    float tfar = 1.0f;
    int hit = 0;
    
    // z position of view plane is -1.0 to fit the cube to the screen quad when axes are aligned, 
    // zoom is -1 and the data set is uniform in each dimension
    // (with FoV of 90° and near plane in range [-1,+1]).
    float4 nearPlanePos = fast_normalize((float4)(imgCoords, -2.0f, 0.0f));
    // transform nearPlane to world space    
    rayDir.x = dot((float4)(viewMat[ 0], viewMat[ 1], viewMat[ 2], viewMat[ 3]), nearPlanePos);
    rayDir.y = dot((float4)(viewMat[ 4], viewMat[ 5], viewMat[ 6], viewMat[ 7]), nearPlanePos);
    rayDir.z = dot((float4)(viewMat[ 8], viewMat[ 9], viewMat[10], viewMat[11]), nearPlanePos);
    rayDir.w = dot((float4)(viewMat[12], viewMat[13], viewMat[14], viewMat[15]), nearPlanePos);
    
    // origin is translation vector of view matrix
    float4 camPos = (float4)(viewMat[3], viewMat[7], viewMat[11], 1.0f);
    rayDir = fast_normalize(rayDir);

    if (orthoCam)
    {
        camPos = (float4)(viewMat[3], viewMat[7], viewMat[11], 1.0f);
        float4 viewPlane_x = (float4)(viewMat[0], viewMat[4], viewMat[ 8], viewMat[12]);
        float4 viewPlane_y = (float4)(viewMat[1], viewMat[5], viewMat[ 9], viewMat[13]);
        float4 viewPlane_z = (float4)(viewMat[2], viewMat[6], viewMat[10], viewMat[14]);
        rayDir = -viewPlane_z;
        rayDir.w = 0;
        rayDir = fast_normalize(rayDir);
        nearPlanePos = camPos + imgCoords.x*viewPlane_x + imgCoords.y*viewPlane_y;
        nearPlanePos *= length(camPos); // TODO fix volume scaling issues
        camPos = nearPlanePos;
    }

    hit = intersectBox(camPos, rayDir, &tnear, &tfar);
    if (!hit)
    {
        write_imagef(outData, texCoords, background);
        return;
    }
    tnear = max(0.0f, tnear);     // clamp to near plane

    float4 result = background;
    float alpha = 0.0f;
    float4 pos = (float4)(0);
    float density = 0.0f;
    float4 tfColor = (float4)(0);
    float opacity = 0.0f;
    float t = 0.0f;

    float3 voxLen = (float3)(2.f) / convert_float3(get_image_dim(volData).xyz);
    float3 entryPoint = camPos.xyz + tnear*rayDir.xyz;
    float3 exitPoint = camPos.xyz + tfar*rayDir.xyz;
    float3 bboxMin = (float3)(-1.f);
    float3 bboxMax = (float3)( 1.f);

    if (checkEdges(entryPoint, bboxMin, bboxMax))
    {
        result.xyz = (float3)(0.f);
        opacity = 1.f;
    }

    uint i = 0;
    // raycasting loop: front to back raycasting with early ray termination
    while (t < tfar)
    {
        t = tnear + stepSize*i;     // recalculate t to avoid numerical drift
        pos = camPos + t*rayDir;
        pos = pos * 0.5f + 0.5f;    // normalize to [0,1]

        density = useLinear ? read_imagef(volData,  linearSmp, pos).x :
                              read_imagef(volData, nearestSmp, pos).x;
        tfColor = read_imagef(tffData, linearSmp, density);         // map density to color
        tfColor.xyz *= useIllum ? illumination(volData, pos) : 1.f; // illumination
        tfColor.xyz = background.xyz - tfColor.xyz;

        // Taylor expansion approximation
        opacity = 1.f - native_powr(1.f - tfColor.w, stepSizeFactor);
        result.xyz = result.xyz - tfColor.xyz * opacity * (1.f - alpha);
        alpha = alpha + opacity * (1.f - alpha);

        if (alpha > ERT_THRESHOLD) break;   // early ray termination check
        ++i;
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

