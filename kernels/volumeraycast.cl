#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

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
int intersectBox(float4 rayOrig, float4 rayDir, float2 splitRange, float *tnear, float *tfar)
{
    // compute intersection of ray with all six bbox planes
    float4 invRay = native_divide((float4)(1.0f, 1.0f, 1.0f, 1.0f), rayDir);
    float4 tBot = invRay * ((float4)(-1.0f, -1.0f, splitRange.x, 1.0f) - rayOrig);
    float4 tTop = invRay * ((float4)(1.0f, 1.0f, splitRange.y, 1.0f) - rayOrig);

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


float illumination(image3d_t vol, const float4 pos, const uint4 volRes)
{
    float3 volResf;
    volResf.x = (float)(volRes.x);
    volResf.y = (float)(volRes.y);
    volResf.z = (float)(volRes.z);
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
    float3 l = fast_normalize((float4)(100.0f, 100.0f, 100.0f, 0.0f) - pos).xyz;

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

int calcSplitId(float2 *imgCoords, uint viewSplit, uint showThumbnails)
{
    int viewId = 0;
    float2 ic = *imgCoords;
    ulong2 splitSize = (ulong2)((0.15f*(float)get_global_size(0)),
                                (0.15f*(float)get_global_size(1)));
    ulong2 invThumbSize = (ulong2)(get_global_size(0) - splitSize.x,
                                   get_global_size(1) - splitSize.y);
    int2 o = (int2)(10, 10);
    // magnify viewSplit id and show corner thumbnails for the others
    if (viewSplit == 0)         // top left
    {
        if (ic.x > invThumbSize.x-o.x && ic.y < splitSize.y+o.y && showThumbnails)  // top right
        {
            ic.x -= (float)invThumbSize.x-o.x;
            ic.y -= o.y;
            viewId = 1;
        }
        else if (ic.x < splitSize.x+o.x && ic.y > invThumbSize.y-o.y && showThumbnails)  // bottom left
        {
            ic.x -= o.x;
            ic.y -= (float)invThumbSize.y-o.y;
            viewId = 2;
        }
        else if (ic.x > invThumbSize.x-o.x && ic.y > invThumbSize.y-o.y && showThumbnails)  // bottom right
        {
            ic -= (float2)(invThumbSize.x-o.x, invThumbSize.y-o.y);
            viewId = 3;
        }
        else
        {
            viewId = 0;
            splitSize = (ulong2)(get_global_size(0), get_global_size(1));
        }
    }
    else if (viewSplit == 1)    // top right
    {
        if (ic.x < splitSize.x+o.x && ic.y < splitSize.y+o.y && showThumbnails)  // top left
        {
            ic.x -= o.x;
            ic.y -= o.y;
            viewId = 0;
        }
        else if (ic.x < splitSize.x+o.x && ic.y > invThumbSize.y-o.y && showThumbnails)  // bottom left
        {
            ic.x -= o.x;
            ic.y -= (float)invThumbSize.y-o.y;
            viewId = 2;
        }
        else if (ic.x > invThumbSize.x-o.x && ic.y > invThumbSize.y-o.y && showThumbnails)  // bottom right
        {
            ic -= (float2)(invThumbSize.x-o.x, invThumbSize.y-o.y);
            viewId = 3;
        }
        else
        {
            viewId = 1;
            splitSize = (ulong2)(get_global_size(0), get_global_size(1));
        }
    }
    else if (viewSplit == 2)    // bottom left
    {
        if (ic.x < splitSize.x+o.x && ic.y < splitSize.y+o.y && showThumbnails)  // top left
        {
            ic.x -= o.x;
            ic.y -= o.y;
            viewId = 0;
        }
        else if (ic.x > invThumbSize.x-o.x && ic.y < splitSize.y+o.y && showThumbnails)  // top right
        {
            ic.x -= (float)invThumbSize.x-o.x;
            ic.y -= o.y;
            viewId = 1;
        }
        else if (ic.x > invThumbSize.x-o.x && ic.y > invThumbSize.y-o.y && showThumbnails)  // bottom right
        {
            ic -= (float2)(invThumbSize.x-o.x, invThumbSize.y-o.y);
            viewId = 3;
        }
        else
        {
            viewId = 2;
            splitSize = (ulong2)(get_global_size(0), get_global_size(1));
        }
    }
    else if (viewSplit == 3)    // bottom right
    {
        if (ic.x < splitSize.x+o.x && ic.y < splitSize.y+o.y && showThumbnails)  // top left
        {
            ic.x -= o.x;
            ic.y -= o.y;
            viewId = 0;
        }
        else if (ic.x > invThumbSize.x-o.x && ic.y < splitSize.y+o.y && showThumbnails)  // top right
        {
            ic.x -= (float)invThumbSize.x-o.x;
            ic.y -= o.y;
            viewId = 1;
        }
        else if (ic.x < splitSize.x+o.x && ic.y > invThumbSize.y-o.y && showThumbnails)  // bottom left
        {
            ic.x -= o.x;
            ic.y -= (float)invThumbSize.y-o.y;
            viewId = 2;
        }
        else
        {
            viewId = 3;
            splitSize = (ulong2)(get_global_size(0), get_global_size(1));
        }
    }
    else // split in 4 equal parts
    {
        splitSize.x = get_global_size(0) / 2ul;
        splitSize.y = get_global_size(1) / 2ul;
        if (ic.x > splitSize.x)     // right
        {
            ic.x -= splitSize.x;
            viewId = 1;             // top right
            if (ic.y > splitSize.y) // bottom right
            {
                ic.y -= splitSize.y;
                viewId = 3;
            }
        }
        else                        // left
        {
            viewId = 0;             // top left
            if (ic.y > splitSize.y) // bottom left
            {
                ic.y -= splitSize.y;
                viewId = 2;
            }
        }
    }

    float aspectRatio = native_divide((float)splitSize.y, (float)(splitSize.x));
    aspectRatio = min(aspectRatio, native_divide((float)splitSize.x, (float)(splitSize.y)));
    int maxSize = splitSize.y; //max(splitSize.x, splitSize.y);

    (*imgCoords).x = native_divide((ic.x + 0.5f), (float)(maxSize)) * 2.0f;
    (*imgCoords).y = native_divide((ic.y + 0.5f), (float)(maxSize)) * 2.0f;
    // calculate correct offset based on aspect ratio
    *imgCoords -= splitSize.x > splitSize.y ?
                        (float2)(1.0f/aspectRatio, 1.0f) : (float2)(aspectRatio, 1.0f);

    // flip y-coordinate to point in right direction
    (*imgCoords).y *= -1.0f;
    return viewId;
}


float4 getReorderedPos(const int viewId, const uint4 volRes, const float4 pos,
                       __constant uint *order)
{
    float4 posReordered;
    uint4 minMax = (uint4)(viewId*volRes.x, (viewId+1)*volRes.x - 1,
                           viewId*volRes.y, (viewId+1)*volRes.y - 1);
    uint2 orderId;
    orderId.x = clamp((uint)(pos.x * (volRes.x-1) + viewId*volRes.x), minMax.x, minMax.y);
    orderId.y = clamp((uint)(pos.y * (volRes.y-1) + viewId*volRes.y), minMax.z, minMax.w);
    posReordered.x = native_divide((float)(order[orderId.x]), (float)volRes.x-1.0f);
    posReordered.y = native_divide((float)(order[orderId.y]), (float)volRes.y-1.0f);
    posReordered.zw = pos.zw;
    return posReordered;
}


float readEdgeLife(uint useMipView, float2 pos,
                   __read_only image2d_t edgeLife0, __read_only image2d_t edgeLife1,
                   __read_only image2d_t edgeLife2, __read_only image2d_t edgeLife3,
                   __read_only image2d_t edgeLife4, __read_only image2d_t edgeLife5,
                   __read_only image2d_t edgeLife6, __read_only image2d_t edgeLife7)
{
    if (useMipView == 0)
        return read_imagef(edgeLife0, nearestSmp, pos).x;
    else if (useMipView == 1)
        return read_imagef(edgeLife1, nearestSmp, pos).x;
    else if (useMipView == 2)
        return read_imagef(edgeLife2, nearestSmp, pos).x;
    else if (useMipView == 3)
        return read_imagef(edgeLife3, nearestSmp, pos).x;
    else if (useMipView == 4)
        return read_imagef(edgeLife4, nearestSmp, pos).x;
    else if (useMipView == 5)
        return read_imagef(edgeLife5, nearestSmp, pos).x;
    else if (useMipView == 6)
        return read_imagef(edgeLife6, nearestSmp, pos).x;
    else if (useMipView == 7)
        return read_imagef(edgeLife7, nearestSmp, pos).x;
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


bool approxEq(float x, float y)
{
    float delta = 0.005f;
    return (x < y + delta) && (y < x + delta);
}


bool approxEq2(float2 a, float2 b)
{
    return approxEq(a.x, b.x) && approxEq(a.y, b.y);
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


__kernel void volumeRender(__read_only image3d_t volData,
                           __write_only image2d_t outData,
                           const uint4 volRes,
                           __constant uint *order,
                           const float4 opacScaling,          // 4 components map to 4 views
                           __read_only image1d_array_t tffData,     // constant transfer function values
                           const float8 tffRange,             // 8 components map to 4 views min+max
                           const uint4 tffMetric,             // 4 components map to 4 views
                           const uint4 tffRepeat,             // 4 components map to 4 views
                           const uint4 useMipLayer,           // 4 components map to 4 views
                           const uint scaleTime,
                           __read_only image3d_t mipLayer1,
                           __read_only image3d_t mipLayer2,
                           __read_only image3d_t mipLayer3,
                           __read_only image3d_t mipLayer4,
                           __read_only image3d_t mipLayer5,
                           __read_only image3d_t mipLayer6,
                           __read_only image3d_t mipLayer7,
                           __read_only image2d_t edgeLife0,
                           __read_only image2d_t edgeLife1,
                           __read_only image2d_t edgeLife2,
                           __read_only image2d_t edgeLife3,
                           __read_only image2d_t edgeLife4,
                           __read_only image2d_t edgeLife5,
                           __read_only image2d_t edgeLife6,
                           __read_only image2d_t edgeLife7,
                           const float4 pickedItemIn,
                           const float4 pickCol,               // ========
                           const float stepSizeFactor,
                           __read_only image2d_t overlay,
                           __constant float *viewMat,
                           __constant float *splitMarks,
                           const float splitDist,
                           const uint splitCnt,
                           const float zScale,
                           const uint viewSplit,
                           const float lodFactor,
                           const uint showThumbnails,
                           const float8 splitRange,           // 8 components map to 4 views min+max
                           const uint4 orthoCam,
                           const int2 pickCoord,
                           __global float4 *pickedItemOut,
                           const int2 lensPos,
                           const int2 lensView,
                           const int2 lensSize,
                           const uint4 useIllumination,
                           const uint4 pickHighlight,
                           __read_only image1d_array_t orderedMetrics, // array of float textures: distance, edge cnt, bandwidth, profile, lin arrangement
                           __read_only image1d_array_t unorderedMetrics,
                           const uint4 tffWeightOpac,
                           const int4 slicePlanes,
                           __global float8 *pickedListOut,
                           __global uint *gCount
                           )
{
    float maxRes = (float)max(volRes.x, max(volRes.y, volRes.z));
    float stepSize = native_divide(stepSizeFactor, maxRes);
    stepSize *= 8.0f; // normalization to octile

    uint idX = (uint)get_global_id(0);
    uint idY = (uint)get_global_id(1);
    long gId = idX + idY * get_global_size(0);
    if (idX >= get_global_size(0) || idY >= get_global_size(1))
        return;
    int2 texCoords = (int2)(idX, idY);

    float2 imgCoords = (float2)(idX, idY);
    int viewId = calcSplitId(&imgCoords, viewSplit, showThumbnails);
    bool isPicked = false;
    bool highlightPicking = false;

    // if cube segment is selected, spawn a second view showing only the split segment
    float2 clampSplit = (float2)(getf8(splitRange, viewId*2), getf8(splitRange, viewId*2 + 1));

    // screen space splitting
    if (viewSplit == 4 &&
            (texCoords.x == get_global_size(0)/2 || texCoords.y == get_global_size(1)/2))
    {
        // draw tile border
        write_imagef(outData, texCoords, (float4)(0.2f, 0.2f, 0.2f, 1.0f));
        return;
    }

    float4 rayDir = (float4)(0, 0, 0, 0);
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

    if (getui4(orthoCam, viewId))
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

    hit = intersectBox(camPos, rayDir, clampSplit, &tnear, &tfar);
    if (!hit)
    {
        // paint overlay background
        write_imagef(outData, texCoords, (float4)(1.0f, 1.0f, 1.0f, 0.0f));
        return;
    }
    tnear = max(0.0f, tnear);     // clamp to near plane

    float initialDist = fast_distance(camPos, (float4)(0));
    float dist = initialDist;
    float minImgRes = (float)min(get_global_size(0), get_global_size(1));
    // factor to account for volume to image resoltion ratio
    dist *= (maxRes / minImgRes);
    dist *= lodFactor;

    float4 color = (float4)(1, 1, 1, 0);
    float4 illumColor = (float4)(0);
    float alpha = 0.0f;
    float4 pos = (float4)(0);
    float4 weight = (float4)(0);
    float4 tfColor = (float4)(0);
    float opacity = 0.0f;
    float t = 0.0f;
    uint i = 0;
    float3 hasData;

    int2 gSize = (int2)(get_global_size(0), get_global_size(1));
    // set lens
    float lensDist = distance((float2)(lensPos.x, lensPos.y), (float2)(texCoords.x, texCoords.y));
    if (lensDist < (float)lensSize.x)
        viewId = lensView.x;
    else if (floor(lensDist) == lensSize.x)  // lens border
    {
        write_imagef(outData, texCoords, (float4)(0.0f, 0.0f, 0.0f, 1.0f));
        return;
    }

    float2 tffRangeView = (float2)(getf8(tffRange, viewId*2), getf8(tffRange, viewId*2 + 1));
    float opacScalingView = getf4(opacScaling, viewId);
    uint repeatView = getui4(tffRepeat, viewId);
    uint tffMetricView = getui4(tffMetric, viewId);
    uint tffWeightOpacView = getui4(tffWeightOpac, viewId);
    uint useMipView = getui4(useMipLayer, viewId);
    uint3 mipRes = (uint3)(volRes.xy/(1<<useMipView), volRes.z/(scaleTime ? (1<<useMipView) : 1));
    float3 pickMarkOffset = (float3)((1 << useMipView)/(float)(volRes.x),
                                     (1 << useMipView)/(float)(volRes.y),
                                     scaleTime ? (1 << useMipView)/(float)(volRes.z) : 1.f/volRes.z);
    uint useIllumView = getui4(useIllumination, viewId);

    // draw bounding box
    pos = camPos + tnear*rayDir;
    pos = pos * 0.5f + 0.5f;
    float3 voxLen = (float3)(2.f / volRes.x, 2.f / volRes.y, 2.f / volRes.z);
    if (checkBoundingBox(pos.xyz, voxLen, (float2)(0, 1)))
    {
        color.xyz = (float3)(0, 0, 0);
        opacity = 1.f;
    }

    float3 entryPoint = camPos.xyz + tnear*rayDir.xyz;
    float3 exitPoint = camPos.xyz + tfar*rayDir.xyz;
    float3 bboxMin = (float3)(-1., -1., clampSplit.x);
    float3 bboxMax = (float3)(1., 1., clampSplit.y);

    if (checkEdges(entryPoint, bboxMin, bboxMax))
    {
        color.xyz = (float3)(0, 0, 0);
        opacity = 1.f;
    }
    float3 posNoReorder;

// MEM ACCESS
//int oldVal = 0;

    // raycasting loop
    // march along ray from front to back, accumulating color
    while (true)
    {
        hasData.x = true;
        hasData.y = true;
        hasData.z = false;
        t = (tnear + stepSize*i);
        //dist = initialDist + t;
        pos = camPos + t*rayDir;
        pos = pos * 0.5f + 0.5f;    // normalize to [0;1]

        posNoReorder = pos.xyz;

        // insert split spaces
        float tmpZ = pos.z*zScale;
        if (clampSplit.x <= -1.f || clampSplit.y >= 1)  // whole volume
        {
            for (int j = 0; j < splitCnt; ++j)
            {
                if (tmpZ > splitMarks[j] + splitDist*j)
                {
                    if (tmpZ < (splitMarks[j] + splitDist*(j)) + voxLen.z ||
                        (tmpZ > (splitMarks[j] + splitDist*(j+1)) - voxLen.z
                        && tmpZ < splitMarks[j] + splitDist*(j+1)))
                        {
                            hasData.y = false;
                            if (approxEq(pos.x, 0) || approxEq(pos.y, 0) ||
                                approxEq(pos.x, 1) || approxEq(pos.y, 1))
                            {
                                color.xyz = (float3)(0, 0, 0);
                                opacity = 1.f;
                            }
                        }
                    if (tmpZ < (splitMarks[j] + splitDist*(j+1)))
                        hasData.x = false;
                    pos.z -= splitDist/zScale;
                }
            }
            pos.z *= zScale;
        }

        // highlight picking after reordering
        uint3 normPos = (uint3)(pos.x * mipRes.x, pos.y * mipRes.y, pos.z * mipRes.z);
        uint3 normPickPos = (uint3)(pickedItemIn.x * mipRes.x, pickedItemIn.y * mipRes.y, pickedItemIn.z * mipRes.z);
        if (   (normPos.x == normPickPos.x && normPos.y == normPickPos.y && pickHighlight.z)
            || (normPos.x == normPickPos.x && normPos.z == normPickPos.z && pickHighlight.y)
            || (normPos.y == normPickPos.y && normPos.z == normPickPos.z && pickHighlight.x)
            || (normPos.x == normPickPos.x && normPos.y == normPickPos.y && normPos.z == normPickPos.z) )
        {
            highlightPicking = true;
            useMipView = getui4(useMipLayer, lensView.y);
        }
        else
        {
            highlightPicking = false;
            useMipView = getui4(useMipLayer, viewId);
        }

        // read_image always returns a 4 component vector: (r, 0.0, 0.0, 1.0)
        // We have to use floating point textures in order to use linear interpolation. (!)
        if (useMipView == 0 || (useMipView > 100 && dist < 2.0f))
        {
            pos = getReorderedPos(viewId, volRes, pos, order);
            weight = read_imagef(volData, nearestSmp, pos).xxxx;
        }
        else if (useMipView == 1)
            weight = read_imagef(mipLayer1, nearestSmp, pos);
        else if (useMipView == 2)
            weight = read_imagef(mipLayer2, nearestSmp, pos);
        else if (useMipView == 3)
            weight = read_imagef(mipLayer3, nearestSmp, pos);
        else if (useMipView == 4)
            weight = read_imagef(mipLayer4, nearestSmp, pos);
        else if (useMipView == 5)
            weight = read_imagef(mipLayer5, nearestSmp, pos);
        else if (useMipView == 6)
            weight = read_imagef(mipLayer6, nearestSmp, pos);
        else if (useMipView == 7)
            weight = read_imagef(mipLayer7, nearestSmp, pos);
        // automatic (distance based) LOD
        else if (useMipView > 100)
        {
            if (dist < 4.0f)
            {
                pos = getReorderedPos(viewId, volRes, pos, order);
                weight = mix(read_imagef(volData, nearestSmp, pos),
                              read_imagef(mipLayer1, nearestSmp, pos), (dist - 2.0f) / 2.0f);
            }
            else if (dist < 8.0f)
                weight = mix(read_imagef(mipLayer1, nearestSmp, pos),
                              read_imagef(mipLayer2, nearestSmp, pos), (dist - 4.0f) / 4.0f);
            else if (dist < 16.0f)
                weight = mix(read_imagef(mipLayer2, nearestSmp, pos),
                              read_imagef(mipLayer3, nearestSmp, pos), (dist - 8.0f) / 8.0f);
            else if (dist < 32.0f)
                weight = mix(read_imagef(mipLayer3, nearestSmp, pos),
                              read_imagef(mipLayer4, nearestSmp, pos), (dist - 16.0f) / 16.0f);
            else if (dist < 64.0f)
                weight = mix(read_imagef(mipLayer4, nearestSmp, pos),
                              read_imagef(mipLayer5, nearestSmp, pos), (dist - 32.0f) / 32.0f);
            else if (dist < 128.0f)
                weight = mix(read_imagef(mipLayer5, nearestSmp, pos),
                              read_imagef(mipLayer6, nearestSmp, pos), (dist - 64.0f) / 64.0f);
            else
                weight = read_imagef(mipLayer6, nearestSmp, pos);
        }
        float viewWeight = getf4(weight, viewId);
        viewWeight = clamp(viewWeight, tffRangeView.x, tffRangeView.y);
        viewWeight *= 1.f/fabs(tffRangeView.y - tffRangeView.x);
        float tffPos = fmod(viewWeight*repeatView, 1);
        // weight based tff
        tfColor.w = read_imagef(tffData, linearSmp, (float2)(tffPos, viewId)).w;

        // apply transfer function
        if (tffMetricView == 1) // time based tff
            tffPos = fmod(pos.z*repeatView, 1);
        else if (tffMetricView == 2) // edge lifetime based tff
        {
            tffPos = readEdgeLife(useMipView, pos.xy, edgeLife0, edgeLife1, edgeLife2, edgeLife3,
                                  edgeLife4, edgeLife5, edgeLife6, edgeLife7);
            tffPos = clamp(tffPos, tffRangeView.x, tffRangeView.y);
            tffPos *= 1.f/fabs(tffRangeView.y - tffRangeView.x);
        }
        else if (tffMetricView == 3) // reordering based tff
        {
            tffPos = order[(int)(pos.x * (volRes.x-1) + viewId*volRes.x)]/(float)(volRes.x - 1);
        }
        else if (tffMetricView == 4)    // distance
            tffPos = read_imagef(unorderedMetrics, nearestSmp, (float2)(pos.z, 0)).x;
        else if (tffMetricView == 5)    // edge count
            tffPos = read_imagef(unorderedMetrics, nearestSmp, (float2)(pos.z, 1)).x;
        else if (tffMetricView == 6)    // bandwidth
            tffPos = read_imagef(orderedMetrics, nearestSmp, (float2)(pos.z, 0)).x;
        else if (tffMetricView == 7)    // profile
            tffPos = read_imagef(orderedMetrics, nearestSmp, (float2)(pos.z, 1)).x;
        else if (tffMetricView == 8)    // linear arrangement
            tffPos = read_imagef(orderedMetrics, nearestSmp, (float2)(pos.z, 2)).x;

        if (tfColor.w)
        {
            tfColor.xyz = read_imagef(tffData, linearSmp, (float2)(tffPos, viewId)).xyz;
            if (!tffWeightOpacView)
                tfColor.w = read_imagef(tffData, linearSmp, (float2)(tffPos, viewId)).w;
        }

        // illumination
        if (useIllumView)
        {
            if (useMipView == 1)
                tfColor.xyz = mix(tfColor.xyz, tfColor.xyz
                                    * illumination(mipLayer1, pos, volRes/(1<<useMipView)), 0.5f);
            else if (useMipView == 2)
                tfColor.xyz = mix(tfColor.xyz, tfColor.xyz
                                    * illumination(mipLayer2, pos, volRes/(1<<useMipView)), 0.5f);
            else if (useMipView == 3)
                tfColor.xyz = mix(tfColor.xyz, tfColor.xyz
                                    * illumination(mipLayer3, pos, volRes/(1<<useMipView)), 0.5f);
            else if (useMipView == 4)
                tfColor.xyz = mix(tfColor.xyz, tfColor.xyz
                                    * illumination(mipLayer4, pos, volRes/(1<<useMipView)), 0.5f);
            else if (useMipView == 5)
                tfColor.xyz = mix(tfColor.xyz, tfColor.xyz
                                    * illumination(mipLayer5, pos, volRes/(1<<useMipView)), 0.5f);
            else if (useMipView == 6)
                tfColor.xyz = mix(tfColor.xyz, tfColor.xyz
                                    * illumination(mipLayer6, pos, volRes/(1<<useMipView)), 0.5f);
        }

        // cull split distance and scale opacity (user input)
        tfColor.w = clamp(tfColor.w * opacScalingView, 0.0f, 1.0f); // space between cuts
        if (!hasData.x)
        {
            tfColor.xyz = (float3)(0.9f, 0.9f, 0.9f);
            tfColor.w = 0.01f;
        }
//        if (!hasData.y)
//        {
//            tfColor.xyz = (float3)(0.8f, 0.8f, 0.8f);
//            tfColor.w = 0.1f;
//        }

        // add to pick list
        float weightOriginal = read_imagef(volData, nearestSmp, pos).x;
        if (highlightPicking && weightOriginal &&
            normPos.x == normPickPos.x &&
            normPos.y == normPickPos.y &&
            normPos.z == normPickPos.z)
        {
            float8 pickedOut = (float8)(pos.xyz, weightOriginal,
                                        read_imagef(edgeLife0, nearestSmp, pos.xy).x, 0,0,0);
            int buffId = floor((pos.x - normPos.x/(float)(mipRes.x))* 127.f);      // x-coord
            buffId += floor((pos.y - normPos.y/(float)(mipRes.y)) * 127.f) * 128;  // y-coord
            if (buffId < 128*128)
                pickedListOut[buffId] = pickedOut;
        }

        // highlight picking
        if (highlightPicking && tfColor.w > 0)
        {
            tfColor.w = clamp(tfColor.w * (1.f + (float)(pickHighlight.w)*0.1f), 0.f, 1.f);
            tfColor.xyz = mix(tfColor.xyz, pickCol.xyz, 0.5f);
        }
        else //if (pickHighlight.x || pickHighlight.y || pickHighlight.z)
            tfColor.w *= (1.f/clamp((float)(pickHighlight.w)*0.1f, 1.0f, 100.f));

        // check for picking (1st raycast pass)
        if ((!isPicked) && (pickCoord.x == texCoords.x) && (pickCoord.y == texCoords.y)
            && (tfColor.w > 0))
        {
            float8 picked = (float8)(pos.xyz, getf4(weight, viewId),
                                     read_imagef(edgeLife0, nearestSmp, pos.xy).x, 0,0,0);
            pickedItemOut[0] = picked.xyzw;
            isPicked = true;
            return;
        }

        if ((slicePlanes.z >= 0) && approxEq(pos.z, (float)(slicePlanes.z)/(float)(volRes.z)))
        {
            tfColor.xyz = (float3)(0.8f, 0.8f, 0.8f);
            tfColor.w = 0.5f;
        }

        opacity = 1.0f - pow(1.0f - tfColor.w, stepSizeFactor);
        color.xyz = color.xyz - ((float3)(1.0f, 1.0f, 1.0f) - tfColor.xyz) * opacity * (1.0f - alpha);
        alpha = alpha + opacity * (1.0f - alpha);

        //oldVal = atomic_inc(gCount); // MEM ACCESS
        if (alpha > 0.98f) break;   // ERT
        if (t > tfar) break;        // out of box
        ++i;
    }

// MEM ACCESS
//write_imagef(outData, texCoords, (float4)((float)(oldVal / (double)(get_global_size(0)*get_global_size(1)*i) / 1.f), 0, 0, 1));
//return;


    if (checkBoundingBox(posNoReorder.xyz, voxLen, (float2)(0, 1)))
    {
        color.xyz = (float3)(0, 0, 0);
        opacity = 1.f;
    }
    if (checkEdges(exitPoint, bboxMin, bboxMax))
    {
        color.xyz = (float3)(0, 0, 0);
        opacity = 1.f;
    }
    // TODO add overlay color
    if (!(clampSplit.x > -1.f || clampSplit.y < 1) && !(showThumbnails && (viewId != viewSplit)))
    {
        float4 overlayCol = read_imagef(overlay, nearestSmp, texCoords);
        color.xyz = color.xyz + (-overlayCol.xyz)*(overlayCol.w) * (1-alpha);
        color.w = alpha;
    }

    write_imagef(outData, texCoords, color);
}


/**
 * Render a volume slice
 */
__kernel void sliceRender(__read_only image3d_t volData,
                          __write_only image2d_t outData,
                          const uint4 volRes,
                          __constant uint *order,
                          const float opacityScaling,
                          __read_only image1d_array_t tffData,
                          const float8 tffRange,
                          const uint4 tffMetric,             // 4 components map to 4 views
                          const uint4 tffRepeat,             // 4 components map to 4 views
                          const uint4 useMipLayer,
                          const uint scaleTime,
                          __read_only image3d_t mipLayer1,
                          __read_only image3d_t mipLayer2,
                          __read_only image3d_t mipLayer3,
                          __read_only image3d_t mipLayer4,
                          __read_only image3d_t mipLayer5,
                          __read_only image3d_t mipLayer6,
                          __read_only image3d_t mipLayer7,
                          __read_only image2d_t edgeLife0,
                          __read_only image2d_t edgeLife1,
                          __read_only image2d_t edgeLife2,
                          __read_only image2d_t edgeLife3,
                          __read_only image2d_t edgeLife4,
                          __read_only image2d_t edgeLife5,
                          __read_only image2d_t edgeLife6,
                          __read_only image2d_t edgeLife7,
                          const float4 pickedItemIn,
                          const float4 pickCol,
                          const uint axis,
                          const uint sliceId,
                          const uint viewId
//                          const uint scaling,
                          )
{
    uint useMipView = getui4(useMipLayer, viewId);
    uint idX = get_global_id(0);
    uint idY = get_global_id(1);
    float2 tffRangeView = (float2)(getf8(tffRange, viewId*2), getf8(tffRange, viewId*2 + 1));
    uint repeatView = getui4(tffRepeat, viewId);
    uint tffMetricView = getui4(tffMetric, viewId);

    int2 texCoords = (int2)(get_global_id(0), (get_global_size(1)-1) - get_global_id(1));
    float4 tfColor = (float4)(1.0f, 1.0f, 1.0f, 0.0f);

    int4 coord = (int4)(0);
    if (axis == 0)   // x-axis
        coord = (int4)(sliceId/(1<<useMipView), idY/(1<<useMipView), idX/(scaleTime ? (1<<useMipView) : 1), 0);
    if (axis == 1)   // y-axis
        coord = (int4)(idY/(1<<useMipView), sliceId/(1<<useMipView), idX/(scaleTime ? (1<<useMipView) : 1), 0);
    if (axis == 2)   // z-axis
        coord = (int4)(idX, idY, sliceId/(scaleTime ? (1<<useMipView) : 1), 0);

    float4 weight;
    if (useMipView == 0)
    {
        uint2 orderId = (uint2)(idX + viewId*volRes.x, idY + viewId*volRes.y);
        orderId = clamp(orderId, (uint2)(0, 0), (uint2)((viewId+1)*volRes.x-1, (viewId+1)*volRes.y-1));
        if (axis == 0)   // x-axis
            coord = (int4)(order[sliceId], order[orderId.y], idX, 0);
        if (axis == 1)   // y-axis
            coord = (int4)(order[orderId.y], order[sliceId], idX, 0);
        if (axis == 2)   // z-axis
            coord = (int4)(order[orderId.x], order[orderId.y], sliceId, 0);
        weight = read_imagef(volData, nearestSmp, coord).xxxx;
    }
    else if (useMipView == 1)
        weight = read_imagef(mipLayer1, nearestSmp, coord);
    else if (useMipView == 2)
        weight = read_imagef(mipLayer2, nearestSmp, coord);
    else if (useMipView == 3)
        weight = read_imagef(mipLayer3, nearestSmp, coord);
    else if (useMipView == 4)
        weight = read_imagef(mipLayer4, nearestSmp, coord);
    else if (useMipView == 5)
        weight = read_imagef(mipLayer5, nearestSmp, coord);
    else if (useMipView == 6)
        weight = read_imagef(mipLayer6, nearestSmp, coord);
    else if (useMipView == 7)
        weight = read_imagef(mipLayer7, nearestSmp, coord);

    float4 pos = (float4)(coord.x/(float)(volRes.x/(1<<useMipView)),
                          coord.y/(float)(volRes.y/(1<<useMipView)),
                          coord.z/(float)(volRes.z),
                          coord.w);

    float viewWeight = getf4(weight, viewId);
    viewWeight = clamp(viewWeight, tffRangeView.x, tffRangeView.y);
    viewWeight *= 1.f/fabs(tffRangeView.y - tffRangeView.x);
    // weight based tff
    float tffPos = fmod(viewWeight*repeatView, 1);
    tfColor.w = read_imagef(tffData, linearSmp, (float2)(tffPos, viewId)).w;

    // apply transfer function
    if (tffMetricView == 1) // time based tff
        tffPos = fmod(pos.z*repeatView, 1);
    else if (tffMetricView == 2) // edge lifetime based tff
    {
        tffPos = readEdgeLife(useMipView, pos.xy, edgeLife0, edgeLife1, edgeLife2, edgeLife3,
                              edgeLife4, edgeLife5, edgeLife6, edgeLife7);
        tffPos = clamp(viewWeight, tffRangeView.x, tffRangeView.y);
        tffPos *= 1.f/fabs(tffRangeView.y - tffRangeView.x);
    }
    // render white if fully transparent
    if (tfColor.w > 0)
    {
        tfColor.xyz = read_imagef(tffData, linearSmp, (float2)(tffPos, viewId)).xyz;
        tfColor.w *= 10.f;
        tfColor.w *= opacityScaling;
    }
    else
        tfColor.w = 1;
    write_imagef(outData, texCoords, tfColor);
}


/**
 * Scale volume data
 */
__kernel void scaleVolume(__read_only image3d_t inVol,
                          __write_only image3d_t outVol,
                          const uint4 inRes,
                          __constant uint *order,
                          const float4 opacScaling,
                          const uint3 brickSize,
                          const int useScaling,
                          const int fromOriginalData,
                          const uint minMetric,
                          const uint maxMetric,
                          const uint avgMetric,
                          const uint densityMetric)
{
    uint3 coords = (uint3)(get_global_id(0), get_global_id(1), get_global_id(2));
    uint3 brickStart = brickSize * coords;
    uint3 brickEnd = min(brickStart + brickSize - (uint3)(1), inRes.xyz - (uint3)(1));

    float4 densities = (float4)(0);
    float4 minDens = (float4)(1);
    float4 maxDens = (float4)(0);
    float4 avgDens = (float4)(0);
    int cnt = 0;
    float4 edgeDensity = (float4)(0);

    // process brick
    for (uint k = brickStart.z; k <= brickEnd.z; ++k)
    {
        for (uint j = brickStart.y; j <= brickEnd.y; ++j)
        {
            for (uint i = brickStart.x; i <= brickEnd.x; ++i)
            {
                if (fromOriginalData)
                {
                    // 4 orderings on the initial data and pass it to 4-component mimaps
                    uint offset = inRes.x;
                    densities.x = read_imagef(inVol, nearestSmp,
                                    (int4)(order[i], order[j], k, 0)).x;
                    densities.y = read_imagef(inVol, nearestSmp,
                                    (int4)(order[i + offset*1], order[j + offset*1], k, 0)).x;
                    densities.z = read_imagef(inVol, nearestSmp,
                                    (int4)(order[i + offset*2], order[j + offset*2], k, 0)).x;
                    densities.w = read_imagef(inVol, nearestSmp,
                                    (int4)(order[i + offset*3], order[j + offset*3], k, 0)).x;
                }
                else
                    densities = read_imagef(inVol, nearestSmp, (int4)(i, j, k, 0));

                edgeDensity += densities > 0.f ? (float4)(1) : (float4)(0);
                minDens = min(minDens, densities);
                maxDens = max(maxDens, densities);
                avgDens += densities;
                cnt++;
            }
        }
    }

    avgDens /= (float4)(cnt);
    edgeDensity /= (float4)(cnt);

    float4 metric = (float4)(0);
    if (minMetric)
        metric = minDens;
    else if (maxMetric)
        metric = maxDens;
    else if (avgMetric)
        metric = avgDens;
    else if (densityMetric)
        metric = edgeDensity;
    if (useScaling)
        metric *= opacScaling;

    write_imagef(outVol, (int4)(coords.x, coords.y, coords.z, 0), metric);
}



/**
 * Calculate graph metrics
 */
__kernel void calcGraphMetrics(__read_only image3d_t inVol,
                               __global uint *densityBuf,
                               __global uint *edgeCntBuf,
                               const uint4 res,
                               __global uint *bandwidthBuf,
                               __global uint *profileBuf,
                               __global uint *linArrangementBuf,
                               __global uint *maxWeightBuf,
                               __constant uint *order)
{
    uint zId = get_global_id(0);

    uint density = 0;
    uint edgeCnt = 0;
    uint bandwidth = 0;
    uint profile = 0;
    uint linArrangement = 0;
    uint maxWeight = 0;

    for (int j = 0; j < res.y; ++j)
    {
        for (int i = 0; j < res.x; ++i)
        {
            uint weight = read_imageui(inVol, nearestSmp, (int4)(order[i], j, zId, 0)).x;
            density += weight;
            if (weight)
            {
                maxWeight = max(weight, maxWeight);
                edgeCnt++;
                bandwidth = max(order[i], bandwidth);
                profile = max(order[i] - j, profile);
                linArrangement = order[i];
            }
        }
    }

    densityBuf[zId] = density;
    edgeCntBuf[zId] = edgeCnt;
    bandwidthBuf[zId] = bandwidth;
    profileBuf[zId] = profile;
    linArrangementBuf[zId] = linArrangement;
    maxWeightBuf[zId] = maxWeight;
}


/**
 *
 */
__kernel void calcEdgeLifetime(__read_only image3d_t inVol,
                               __write_only image2d_t outEdgeLen,
                               const uint4 res)
{
    int2 texCoords = (int2)(get_global_id(0), (get_global_size(1)-1) - get_global_id(1));

    int edgeLifetime = 0;
    for (int i = 0; i < res.z; ++ i)
    {
        float weight = read_imagef(inVol, nearestSmp, (int4)(texCoords.x, texCoords.y, i, 0)).x;
        if (weight)
            edgeLifetime++;
    }

    write_imagef(outEdgeLen, texCoords, edgeLifetime/(float)(res.z));
}
