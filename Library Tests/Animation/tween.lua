--------------------------------------------------------------------------------
-- Tweening Library
--
-- @author   Mayeko
-- @file     tween.lua
-- @version  1.0
--------------------------------------------------------------------------------

-- TODO : Try forward calculation (and check whether it freezes or not)

do

--------------------------------------------------------------------------------
--< Private data >--------------------------------------------------------------
--------------------------------------------------------------------------------
        
        local render = Game:GetService("RunService").RenderStepped;

        local sqrt = math.sqrt;
        local acos = math.acos;
        local asin = math.asin
        local cos = math.cos;
        local sin = math.sin;
        local rad = math.rad;
        local pow = math.pow
        local abs = math.abs
        local cfn = CFrame.new;
        local cfa = CFrame.Angles;
        local vec = Vector3.new;
        local pi = math.pi
        
        local remove_index = table.remove;
        local co_resume = coroutine.resume;
        local co_create = coroutine.create;
        
        local next   = next;
        local pcall  = pcall;
        local unpack = unpack;
        local error  = error;
        local assert = assert;

        local function fassert(condition, err)
                if (not condition) then
                        if (err) then
                                assert(false, err);
                        else 
                                assert(false);
                        end
                end
        end

-- Penner's easing equations >-------------------------------------------------
        -- Source: https://github.com/EmmanuelOga/easing
        -- Curves: http://hosted.zeh.com.br/tweener/docs/en-us/misc/transitions.html

        local function linear(t, b, c, d)
                return c * t / d + b
        end

        local function inQuad(t, b, c, d)
                t = t / d
                return c * pow(t, 2) + b
        end

        local function outQuad(t, b, c, d)
                t = t / d
                return -c * t * (t - 2) + b
        end

        local function inOutQuad(t, b, c, d)
                t = t / d * 2
                if t < 1 then
                        return c / 2 * pow(t, 2) + b
                else
                        return -c / 2 * ((t - 1) * (t - 3) - 1) + b
                end
        end

        local function outInQuad(t, b, c, d)
                if t < d / 2 then
                        return outQuad (t * 2, b, c / 2, d)
                else
                        return inQuad((t * 2) - d, b + c / 2, c / 2, d)
                end
        end

        local function inCubic (t, b, c, d)
                t = t / d
                return c * pow(t, 3) + b
        end

        local function outCubic(t, b, c, d)
                t = t / d - 1
                return c * (pow(t, 3) + 1) + b
        end

        local function inOutCubic(t, b, c, d)
                t = t / d * 2
                if t < 1 then
                        return c / 2 * t * t * t + b
                else
                        t = t - 2
                        return c / 2 * (t * t * t + 2) + b
                end
        end

        local function outInCubic(t, b, c, d)
                if t < d / 2 then
                        return outCubic(t * 2, b, c / 2, d)
                else
                        return inCubic((t * 2) - d, b + c / 2, c / 2, d)
                end
        end

        local function inQuart(t, b, c, d)
                t = t / d
                return c * pow(t, 4) + b
        end

        local function outQuart(t, b, c, d)
                t = t / d - 1
                return -c * (pow(t, 4) - 1) + b
        end

        local function inOutQuart(t, b, c, d)
                t = t / d * 2
                if t < 1 then
                        return c / 2 * pow(t, 4) + b
                else
                        t = t - 2
                        return -c / 2 * (pow(t, 4) - 2) + b
                end
        end

        local function outInQuart(t, b, c, d)
                if t < d / 2 then
                        return outQuart(t * 2, b, c / 2, d)
                else
                        return inQuart((t * 2) - d, b + c / 2, c / 2, d)
                end
        end

        local function inQuint(t, b, c, d)
                t = t / d
                return c * pow(t, 5) + b
        end

        local function outQuint(t, b, c, d)
                t = t / d - 1
                return c * (pow(t, 5) + 1) + b
        end

        local function inOutQuint(t, b, c, d)
                t = t / d * 2
                if t < 1 then
                        return c / 2 * pow(t, 5) + b
                else
                        t = t - 2
                        return c / 2 * (pow(t, 5) + 2) + b
                end
        end

        local function outInQuint(t, b, c, d)
                if t < d / 2 then
                        return outQuint(t * 2, b, c / 2, d)
                else
                        return inQuint((t * 2) - d, b + c / 2, c / 2, d)
                end
        end

        local function inSine(t, b, c, d)
                return -c * cos(t / d * (pi / 2)) + c + b
        end

        local function outSine(t, b, c, d)
                return c * sin(t / d * (pi / 2)) + b
        end

        local function inOutSine(t, b, c, d)
                return -c / 2 * (cos(pi * t / d) - 1) + b
        end

        local function outInSine(t, b, c, d)
                if t < d / 2 then
                        return outSine(t * 2, b, c / 2, d)
                else
                        return inSine((t * 2) -d, b + c / 2, c / 2, d)
                end
        end

        local function inExpo(t, b, c, d)
                if t == 0 then
                        return b
                else
                        return c * pow(2, 10 * (t / d - 1)) + b - c * 0.001
                end
        end

        local function outExpo(t, b, c, d)
                if t == d then
                        return b + c
                else
                        return c * 1.001 * (-pow(2, -10 * t / d) + 1) + b
                end
        end

        local function inOutExpo(t, b, c, d)
                if t == 0 then return b end
                if t == d then return b + c end
                t = t / d * 2
                if t < 1 then
                        return c / 2 * pow(2, 10 * (t - 1)) + b - c * 0.0005
                else
                        t = t - 1
                        return c / 2 * 1.0005 * (-pow(2, -10 * t) + 2) + b
                end
        end

        local function outInExpo(t, b, c, d)
                if t < d / 2 then
                        return outExpo(t * 2, b, c / 2, d)
                else
                        return inExpo((t * 2) - d, b + c / 2, c / 2, d)
                end
        end

        local function inCirc(t, b, c, d)
                t = t / d
                return(-c * (sqrt(1 - pow(t, 2)) - 1) + b)
        end

        local function outCirc(t, b, c, d)
                t = t / d - 1
                return(c * sqrt(1 - pow(t, 2)) + b)
        end

        local function inOutCirc(t, b, c, d)
                t = t / d * 2
                if t < 1 then
                        return -c / 2 * (sqrt(1 - t * t) - 1) + b
                else
                        t = t - 2
                        return c / 2 * (sqrt(1 - t * t) + 1) + b
                end
        end

        local function outInCirc(t, b, c, d)
                if t < d / 2 then
                        return outCirc(t * 2, b, c / 2, d)
                else
                        return inCirc((t * 2) - d, b + c / 2, c / 2, d)
                end
        end

        local function inElastic(t, b, c, d, a, p)
                if t == 0 then return b end

                t = t / d

                if t == 1        then return b + c end

                if not p then p = d * 0.3 end

                local s

                if not a or a < abs(c) then
                        a = c
                        s = p / 4
                else
                        s = p / (2 * pi) * asin(c/a)
                end

                t = t - 1

                return -(a * pow(2, 10 * t) * sin((t * d - s) * (2 * pi) / p)) + b
        end

        -- a: amplitud
        -- p: period
        local function outElastic(t, b, c, d, a, p)
                if t == 0 then return b end

                t = t / d

                if t == 1 then return b + c end

                if not p then p = d * 0.3 end

                local s

                if not a or a < abs(c) then
                        a = c
                        s = p / 4
                else
                        s = p / (2 * pi) * asin(c/a)
                end

                return a * pow(2, -10 * t) * sin((t * d - s) * (2 * pi) / p) + c + b
        end

        -- p = period
        -- a = amplitud
        local function inOutElastic(t, b, c, d, a, p)
                if t == 0 then return b end

                t = t / d * 2

                if t == 2 then return b + c end

                if not p then p = d * (0.3 * 1.5) end
                if not a then a = 0 end

                local s

                if not a or a < abs(c) then
                        a = c
                        s = p / 4
                else
                        s = p / (2 * pi) * asin(c / a)
                end

                if t < 1 then
                        t = t - 1
                        return -0.5 * (a * pow(2, 10 * t) * sin((t * d - s) * (2 * pi) / p)) + b
                else
                        t = t - 1
                        return a * pow(2, -10 * t) * sin((t * d - s) * (2 * pi) / p ) * 0.5 + c + b
                end
        end

        -- a: amplitud
        -- p: period
        local function outInElastic(t, b, c, d, a, p)
                if t < d / 2 then
                        return outElastic(t * 2, b, c / 2, d, a, p)
                else
                        return inElastic((t * 2) - d, b + c / 2, c / 2, d, a, p)
                end
        end

        local function inBack(t, b, c, d, s)
                if not s then s = 1.70158 end
                t = t / d
                return c * t * t * ((s + 1) * t - s) + b
        end

        local function outBack(t, b, c, d, s)
                if not s then s = 1.70158 end
                t = t / d - 1
                return c * (t * t * ((s + 1) * t + s) + 1) + b
        end

        local function inOutBack(t, b, c, d, s)
                if not s then s = 1.70158 end
                s = s * 1.525
                t = t / d * 2
                if t < 1 then
                        return c / 2 * (t * t * ((s + 1) * t - s)) + b
                else
                        t = t - 2
                        return c / 2 * (t * t * ((s + 1) * t + s) + 2) + b
                end
        end

        local function outInBack(t, b, c, d, s)
                if t < d / 2 then
                        return outBack(t * 2, b, c / 2, d, s)
                else
                        return inBack((t * 2) - d, b + c / 2, c / 2, d, s)
                end
        end

        local function outBounce(t, b, c, d)
                t = t / d
                if t < 1 / 2.75 then
                        return c * (7.5625 * t * t) + b
                elseif t < 2 / 2.75 then
                        t = t - (1.5 / 2.75)
                        return c * (7.5625 * t * t + 0.75) + b
                elseif t < 2.5 / 2.75 then
                        t = t - (2.25 / 2.75)
                        return c * (7.5625 * t * t + 0.9375) + b
                else
                        t = t - (2.625 / 2.75)
                        return c * (7.5625 * t * t + 0.984375) + b
                end
        end

        local function inBounce(t, b, c, d)
                return c - outBounce(d - t, 0, c, d) + b
        end

        local function inOutBounce(t, b, c, d)
                if t < d / 2 then
                        return inBounce(t * 2, 0, c, d) * 0.5 + b
                else
                        return outBounce(t * 2 - d, 0, c, d) * 0.5 + c * .5 + b
                end
        end

        local function outInBounce(t, b, c, d)
                if t < d / 2 then
                        return outBounce(t * 2, b, c / 2, d)
                else
                        return inBounce((t * 2) - d, b + c / 2, c / 2, d)
                end
        end

        -- CFrame.Angles with angles in degrees

        local function deg(x, y, z)
                return cfa(rad(x), rad(y), rad(z));
        end


        -- Convert euler angles to quaternion
        -- Source: [ADAPTED CODE] http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/

        local function quat(heading, attitude, bank)
                if (not heading and not attitude and not bank) then
                        return {0, 0, 0, 1};
                end

                local heading  = rad(heading);
                local attitude = rad(attitude);
                local bank         = rad(bank);

                local c1 = cos(.5*heading);
                local s1 = sin(.5*heading);
                local c2 = cos(.5*attitude);
                local s2 = sin(.5*attitude);
                local c3 = cos(.5*bank);
                local s3 = sin(.5*bank);

                local c1c2 = c1*c2;
                local s1s2 = s1*s2;

                return {
                        c1c2*s3  + s1s2*c3,
                        s1*c2*c3 + c1*s2*s3,
                        c1*s2*c3 - s1*c2*s3,
                        c1c2*c3  - s1s2*s3
                };
        end


        -- Get quaternion from CFrame
        -- Source: [ADAPTED CODE] http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
        local function GetQuaternion(CFrame)
                local px,  py,  pz,
                          m00, m01, m02,
                          m10, m11, m12,
                          m20, m21, m22 = CFrame:components();

                local tr = m00 + m11 + m22;
                local s, qw, qx, qy, qz;

                if tr > 0 then
                        s = sqrt(tr + 1) * 2; -- s = 4*qw 
                        qw = .25 * s;
                        qx = (m21 - m12) / s;
                        qy = (m02 - m20) / s;
                        qz = (m10 - m01) / s;
                elseif (m00 > m11) and (m00 > m22) then 
                        s = sqrt(1 + m00 - m11 - m22) * 2; -- s = 4*qx 
                        qw = (m21 - m12) / s;
                        qx = .25 * s;
                        qy = (m01 + m10) / s;
                        qz = (m02 + m20) / s;
                elseif (m11 > m22) then 
                        s = sqrt(1 + m11 - m00 - m22) * 2; -- s = 4*qy
                        qw = (m02 - m20) / s;
                        qx = (m01 + m10) / s;
                        qy = .25 * s;
                        qz = (m12 + m21) / s;
                else
                        s = sqrt(1 + m22 - m00 - m11) * 2; -- s = 4*qz
                        qw = (m10 - m01) / s;
                        qx = (m02 + m20) / s;
                        qy = (m12 + m21) / s;
                        qz = .25 * s;
                end

                return qx, qy, qz, qw;
        end

        local function MulQuaternion(ax, ay, az, aw, bx, by, bz, bw)
                return  ax * bw + ay * bz - az * by + aw * bx,
                           -ax * bz + ay * bw + az * bx + aw * by,
                                ax * by - ay * bx + az * bw + aw * bz,
                           -ax * bx - ay * by - az * bz + aw * bw;
        end

        -- Slerp two quaternions
        -- Source: [ADAPTED CODE] http://www.geometrictools.com/LibMathematics/Algebra/Wm5Quaternion.inl

        local function SlerpQuaternion(ax, ay, az, aw, bx, by, bz, bw, t) 
                local dot = ax*bx + ay*by + az*bz + aw*bw;

                if (abs(dot) >= 1) then
                        return ax, ay, az, aw;
                end 

                local theta = acos(dot);

                local invSin = 1 / sin(theta);
                local tAngle = t*theta;
                local coeff0 = invSin*sin(theta - tAngle);
                local coeff1 = invSin*sin(tAngle);

                return coeff0*ax + coeff1*bx,
                       coeff0*ay + coeff1*by,
                       coeff0*az + coeff1*bz,
                       coeff0*aw + coeff1*bw;
        end

--------------------------------------------------------------------------------
--< Main >----------------------------------------------------------------------
--------------------------------------------------------------------------------

        local slerp_routines, slerp_waiting = {}, {};
        local lerp_routines,  lerp_waiting  = {}, {};


        -- Runs a new thread or create a pending thread
        -- Private function: called undirectly through _G.tween.slerp.new() -> slerp_new_thread()
        --                                          or _G.tween.lerp.new()  -> lerp_new_thread()

        local function new_thread(callback, routines, waiting, ...)
                if (routines[...]) then
                        waiting[#waiting+1] = {...}; -- adds into pending list
                else -- otherwise runs a new thread
                        local _, assertion, err = pcall(co_resume, co_create(callback), ...);
                        fassert(assertion, err);
                end
        end

        -- Slerp part from "cf_start" to "vec_stop" and "quat_stop" in "frames" frames
        -- Private function: called undirectly through _G.tween.slerp.new() -> slerp_new_thread() -> new_thread()

        local function slerp_new(object, cf_start, vec_stop, quat_stop, frames, easing, p1, p2)
                fassert(object, "attempt to call local 'object' (a nil value)");
                fassert(vec_stop or quat_stop, "attempt to call local 'vec_stop' and 'quat_stop' (nil values)");

                if (slerp_routines[object]) then return nil; end
                slerp_routines[object] = true;

                -- default parameters' value
                frames = frames and 1/abs(frames) or .1; -- default = 10 frames
                easing = easing or  linear;
                if (not cf_start) then
                                if (object:IsA("BasePart"))  then cf_start = object.CFrame;
                        elseif (object:IsA("JointInstance")) then cf_start = object.C0:inverse();
                        end
                end

                -- start data
                local vec_start                  = cf_start.p;
                local aqx, aqy, aqz, aqw = GetQuaternion(cf_start);
                local apx, apy, apz          = vec_start.x, vec_start.y, vec_start.z;

                -- stop data
                local bpx, bpy, bpz;
                local bqx, bqy, bqz, bqw;
                if  (vec_stop) then bpx, bpy, bpz =  vec_stop.x,  vec_stop.y,  vec_stop.z;
                                           else bpx, bpy, bpz = vec_start.x, vec_start.y, vec_start.z;
                end
                if (quat_stop) then bqx, bqy, bqz, bqw = unpack(quat_stop);
                                           else bqx, bqy, bqz, bqw = aqx, aqy, aqz, aqw;
                end

                -- interpolate
                if (object:IsA("BasePart")) then
                        for f = 0, 1, frames do
                                if (not slerp_routines[object]) then break; end
                                local t  = easing(f, 0, 1, 1, p1, p2);
                                local rt = 1-t;
                                object.CFrame = cfn( apx*rt + bpx*t, apy*rt + bpy*t, apz*rt + bpz*t,
                                                                         SlerpQuaternion( aqx, aqy, aqz, aqw,
                                                                                          bqx, bqy, bqz, bqw, t ));
                                render:wait();
                        end
                        if (slerp_routines[object]) then
                                object.CFrame = cfn(bpx, bpy, bpz, bqx, bqy, bqz, bqw); -- t == 1;
                        end
                elseif (object:IsA("JointInstance")) then
                        for f = 0, 1, frames do
                                if (not slerp_routines[object]) then break; end
                                local t  = easing(f, 0, 1, 1, p1, p2);
                                local rt = 1-t;
                                object.C0 = cfn( apx*rt + bpx*t, apy*rt + bpy*t, apz*rt + bpz*t,
                                                                 SlerpQuaternion( aqx, aqy, aqz, aqw,
                                                                                  bqx, bqy, bqz, bqw, t )):inverse();
                                render:wait();
                        end
                        if (slerp_routines[object]) then
                                object.C0 = cfn(bpx, bpy, bpz, bqx, bqy, bqz, bqw):inverse(); -- t == 1;
                        end
                end

                if (slerp_routines[object]) then -- if slerp wasn't stopped
                        for i = 1, #slerp_waiting do -- checks for any waiting thread...
                                local t_waiting = slerp_waiting[i]
                                if (t_waiting[1] == object) then -- ...which will slerp 'object'
                                        slerp_routines[object] = nil;
                                        local _, assertion, err = pcall(co_resume, co_create(slerp_new), unpack(t_waiting, 1, 8)); -- 9 arguments (avoid table.maxn)
                                        fassert(assertion, err);
                                        remove_index(slerp_waiting, i); -- shifts other pending threads
                                        break;
                                end
                        end
                end

        end


        -- Generic lerp
        -- Private function: called undirectly through _G.tween.lerp.new() -> lerp_new_thread() -> new_thread()

        -- Index checker for lerp_new
        local function pcall_index(object, parameter)
                object[parameter] = object[parameter];
        end

        local function lerp_new(object, parameter, callback, frames, easing, p1, p2)
                fassert(object,        "attempt to call local 'object' (a nil value)");
                fassert(parameter, "attempt to call local 'parameter' (a nil value)");
                fassert(callback,  "attempt to call local 'callback' (a nil value)");

                if (lerp_routines[object] or not pcall(pcall_index, object, parameter)) then return nil; end
                lerp_routines[object] = true;
                local a = 3;
                -- default parameters' value
                frames = frames and 1/abs(frames) or .1; -- default = 10 frames
                easing = easing or  linear;
                start  = start  or  object[parameter];

                for f = 0, 1, frames do
                        if (not lerp_routines[object]) then break; end
                        object[parameter] = callback(easing(f, 0, 1, 1, p1, p2), object[parameter]);
                        render:wait();
                end
                object[parameter] = callback(1, object[parameter]);
                
                if (lerp_routines[object]) then -- if lerp wasn't stoppedh
                        for i = 1, #lerp_waiting do -- checks for any waiting thread...
                                local t_waiting = lerp_waiting[i]
                                if (t_waiting[1] == object) then -- ...which will slerp 'object'
                                        lerp_routines[object] = nil;
                                        local _, assertion, err = pcall(co_resume, co_create(lerp_new), unpack(t_waiting, 1, 7)); -- 7 arguments (avoid table.maxn)
                                        fassert(assertion, err);
                                        remove_index(lerp_waiting, i); -- shifts other pending threads
                                        break;
                                end
                        end
                end

        end


        -- Runs a new thread or create a pending thread

        local function slerp_new_thread(...)
                new_thread(slerp_new, slerp_routines, slerp_waiting, ...);
        end

        local function lerp_new_thread(...)
                new_thread(lerp_new, lerp_routines, lerp_waiting, ...);
        end


        -- Checks for any ongoing thread

        local function running(object, routines, waiting)
                if (object) then
                        return routines[object];
                else
                        return next(routines) and true;
                end
        end

        local function slerp_running(object)
                return running(object, slerp_routines, slerp_waiting);
        end

        local function lerp_running(object)
                return running(object, lerp_routines, lerp_waiting);
        end


        -- Stops and removes any ongoing and pending thread

        local function stop(object, routines, waiting)
                if (object) then
                        routines[object] = nil; -- stop any object specific ongoing thread
                        for i = 1, #waiting do -- delete every object specific pending thread
                                if (waiting[i][1] == object) then
                                        remove_index(waiting, i); -- shifts other pending threads
                                end
                        end
                else
                        for i, v in next, routines do
                                routines[i] = nil; -- stop every routines
                        end
                        for i = 1, #waiting do
                                waiting[i] = nil; -- delete every pending thread
                        end
                end
                wait(.05) -- needed to let any ongoing slerp thread detect the stop signal
        end

        local function slerp_stop(object)
                stop(object, slerp_routines, slerp_waiting);
        end

        local function lerp_stop(object)
                stop(object, lerp_routines, lerp_waiting);
        end



--------------------------------------------------------------------------------
--< Public access >-------------------------------------------------------------
--------------------------------------------------------------------------------

        _G.tween = {
                deg      = deg,
                quat_new = quat,
                get_quat = GetQuaternion,
                slerp   = {
                        new     = slerp_new_thread,
                        running = slerp_running,
                        stop    = slerp_stop
                },
                lerp    = {
                        new     = lerp_new_thread,
                        running = lerp_running,
                        stop    = lerp_stop
                },
                easing  = {
                        linear       = linear,
                        inQuad       = inQuad,
                        outQuad      = outQuad,
                        inOutQuad    = inOutQuad,
                        outInQuad    = outInQuad,
                        inCubic      = inCubic ,
                        outCubic     = outCubic,
                        inOutCubic   = inOutCubic,
                        outInCubic   = outInCubic,
                        inQuart      = inQuart,
                        outQuart     = outQuart,
                        inOutQuart   = inOutQuart,
                        outInQuart   = outInQuart,
                        inQuint      = inQuint,
                        outQuint     = outQuint,
                        inOutQuint   = inOutQuint,
                        outInQuint   = outInQuint,
                        inSine       = inSine,
                        outSine      = outSine,
                        inOutSine    = inOutSine,
                        outInSine    = outInSine,
                        inExpo       = inExpo,
                        outExpo      = outExpo,
                        inOutExpo    = inOutExpo,
                        outInExpo    = outInExpo,
                        inCirc       = inCirc,
                        outCirc      = outCirc,
                        inOutCirc    = inOutCirc,
                        outInCirc    = outInCirc,
                        inElastic    = inElastic,
                        outElastic   = outElastic,
                        inOutElastic = inOutElastic,
                        outInElastic = outInElastic,
                        inBack       = inBack,
                        outBack      = outBack,
                        inOutBack    = inOutBack,
                        outInBack    = outInBack,
                        inBounce     = inBounce,
                        outBounce    = outBounce,
                        inOutBounce  = inOutBounce,
                        outInBounce  = outInBounce
                }
        };
end

return;
