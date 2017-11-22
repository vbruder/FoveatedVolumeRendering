#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

#define ERT_THRESHOLD 0.98

constant sampler_t linearSmp = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_CLAMP |
                                CLK_FILTER_LINEAR;
constant sampler_t nearestSmp = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_CLAMP |
                                CLK_FILTER_NEAREST;

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
                           )
{
    int2 globalId = (int2)(get_global_id(0), get_global_id(1));
    if(any(globalId >= get_image_dim(outData)))
        return;

    // pseudo random number [0,1] for ray offsets to avoid moire patterns
    float iptr;
    float rand = fract(sin(dot(convert_float2(globalId),
                       (float2)(12.9898f, 78.233f))) * 43758.5453f, &iptr);

    float maxRes = (float)max(get_image_dim(volData).x, get_image_dim(volData).z);
    float stepSize = native_divide(8.f, maxRes*samplingRate); // normalization to octile

    int2 texCoords = globalId;
    float2 imgCoords;
    imgCoords.x = native_divide((globalId.x + 0.5f), (float)(get_global_size(0))) * 2.f - 1.f;
    imgCoords.y = native_divide((globalId.y + 0.5f), (float)(get_global_size(1))) * 2.f - 1.f;
    imgCoords.y *= -1.f;   // flip y coord

    float4 rayDir = (float4)(0.f);
    float tnear = 0.f;
    float tfar = 1.f;
    int hit = 0;
    
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

    hit = intersectBox(camPos, rayDir, &tnear, &tfar);
    if (!hit || tfar < 0)
    {
        write_imagef(outData, texCoords, background);
        return;
    }

    tnear = max(0.f, tnear + rand*stepSize); // clamp to near plane and offset by 'random' distance
    float4 result = background;
    float alpha = 0.f;
    float4 pos = (float4)(0);
    float density = 0.f;
    float4 tfColor = (float4)(0);
    float opacity = 0.f;
    float t = tnear;
    int3 volRes = get_image_dim(volData).xyz;
    float3 voxLen = (float3)(1.f) / convert_float3(volRes);
    int3 bricksRes = get_image_dim(volBrickData).xyz;
    float3 brickLen = (float3)(1.f) / convert_float3(bricksRes);
    float refSamplingInterval = 1.f / samplingRate;
    float precisionDiv = 1.f;
    if (get_image_channel_data_type(volData) == CLK_UNORM_INT16)
        precisionDiv = 8.f;

    // DDA initialization
    float4 posDDA = (float4)(0);
    float3 invRay = 1.f/(rayDir.xyz);
    if (rayDir.x == 0) invRay.x = FLT_MAX;
    if (rayDir.y == 0) invRay.y = FLT_MAX;
    if (rayDir.z == 0) invRay.z = FLT_MAX;

    float3 step = (float3)1 * sign(rayDir.xyz);
    const float3 deltaT = fabs(brickLen*invRay);
    float3 voxIncr = (float3)0;

    float3 bbHit = camPos.xyz + tnear * rayDir.xyz;     // [-1; 1]
    float3 bbHitNorm = bbHit * 0.5f + 0.5f;             // [ 0; 1]
    float3 tv = bbHitNorm * convert_float3(bricksRes.xyz); // [ 0; volRes]

    int3 brickId = convert_int3(floor(bbHitNorm * convert_float3(bricksRes.xyz)));
    // assert range [0; bricksRes - 1]
    brickId = clamp(brickId, 0, convert_int3(bricksRes.xyz - 1));

    // correct initial distances to first inner grid boundary
    tv -= convert_float3(brickId);      // [0; 1] (fractual part)
    tv *= brickLen;                     // [0; brickLen]
    tv = (brickLen - tv) * invRay;

    uint i = 0;
    float3 p_exit = (float3)(0);
    float t_exit = 0;
    float skip = 0;
    // raycasting loop: front to back raycasting with early ray termination
    while (t < tfar)
    {
        density = read_imagef(volBrickData, nearestSmp, (int4)(brickId, 0)).y;

        // increment to next brick
        voxIncr.x = (tv.x <= tv.y) && (tv.x <= tv.z) ? 1 : 0;
        voxIncr.y = (tv.y <= tv.x) && (tv.y <= tv.z) ? 1 : 0;
        voxIncr.z = (tv.z <= tv.x) && (tv.z <= tv.y) ? 1 : 0;
        tv += voxIncr * deltaT;
        brickId += convert_int3(voxIncr * step);    // [0; res-1]

        p_exit = tv;
        p_exit -= convert_float3(brickId);      // [0; 1] (fractual part)
        p_exit *= brickLen;                     // [0; brickLen]
        p_exit = (brickLen - t_exit) * invRay;
        t_exit = t + length(p_exit);    // FIXME
        // check if brick contains values -> normal sampling
        //float alphaMax = read_imagef(tffData, linearSmp, density).w;
        if (texCoords.x == 400 && texCoords.y == 400)
        {
            printf("%#v3f", p_exit);
            printf(" %#f\n", t);
            printf(" %#f\n", t_exit);
            printf(" %#f\n", stepSize);
        }
        if (density > 0.f)
        {
            while (t < t_exit && i < 1512)
            {
                // t += tnear + stepSize*i;     // recalculate t to avoid numerical drift
                t += stepSize;
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

                if (alpha > ERT_THRESHOLD) break;   // early ray termination check
                ++i;
            }
            t = t_exit;
        }
        else
        {
            t = t_exit;
            skip += 0.01f;
        }

        if (any(brickId < 0) || any(brickId > bricksRes.xyz-1))
            break;
    }
//write_imagef(outData, texCoords, (float4)(skip, 0, 0, 1.0f));
//return;

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
