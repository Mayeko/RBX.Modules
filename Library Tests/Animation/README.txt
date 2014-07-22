Tweening library


Features:

        - Provides functions to lerp any object property and slerp BasePart and JointInstance classes.
        - Provides neat animations functions using Penner's easing equations.
        - Possibility to rotate over 180 degrees using the _G.anim.quat_new() function.
        - Possibility to manage threads. Here are some explanations of how threads are and can be managed:
                > Any ongoing (s)lerp operates into a coroutine and is registered into a table.
                > Thus, any (s)lerp call won't yield your script.
                > You cannot run two (s)lerp on the same object at the same time.
                > If you try to do so, the following threads will be added into a pending list.
                > When a (s)lerp finishes on an object, the function checks automatically if there's any pending thread specific to this object.
                > You can stop a thread specific to an object or every threads using _G.anim.stop([Instance object]);
                > Stopping a thread will also delete every pending object specific thread, or every pending threads if object == nil.
                > You can check whether there's any ongoing (s)lerp on any or a specific object using _G.anim.running([Instance object]);

________


Functions:

1.      _G.anim.deg(
                int heading,  -- in degrees
                int attitude, -- in degrees
                int bank      -- in degrees
        ); -- returns CFrame:Euler;

        e.g.    local rotation = _G.anim.deg(0, 0, 0);

        ________


2.      _G.anim.quat_new(
                int heading,  -- in degrees
                int attitude, -- in degrees
                int bank      -- in degrees
        ); -- returns CFrame:Quaternion;

        e.g.    local rotation = _G.anim.quat_new(0, 0, 0);

        ________


3.      _G.anim.get_quat(
                CFrame cframe
        ); -- returns int x, int y, int z, int w;

        e.g.    local x, y, z, w = _G.anim.get_quat(CFrame.new());

        ________


4.      _G.anim.lerp.new(
                Instance           object,    -- object to apply interpolation
                String             parameter, -- object's parameter to interpolate
                LuaFunction(int t) callback,  -- function which returns parameter's value
                int                [frames,]  -- duration in frames (default: 10)
                LuaFunction        [easing,]  -- easing function (default: linear)
                int                [p1,]      -- parameter for elastic or back type tweening
                int                [p2]       -- parameter for elastic tweening
        );

        e.g.    _G.anim.lerp.new(
                        Workspace.Part,
                        "Size", 
                        function(theta)
                                return Vector3.new(1, 1, 1)*(1-theta) -- starting size
                                     + Vector3.new(10, 10, 10)*theta  -- ending   size
                        end,
                        100,
                        _G.anim.easing.outBounce
                );

        ________


5.      _G.anim.lerp.running(
                Instance [object]
        ); -- returns boolean: true whether there's any running lerp on 'object', otherwise false
                               if object == nil, checks for any running thread

        e.g.     _G.anim.lerp.running(Workspace.Part);

        ________


6.      _G.anim.lerp.stop(
                Instance [object]
        ); -- stops every running lerp on 'object' and deletes every pending threads of this object
              if object == nil, stops every running lerp and deletes every pending threads

        e.g.    _G.anim.lerp.stop(Workspace.Part);

        ________


7.      _G.anim.slerp.new(
                Instance           object,    -- object to interpolate
                CFrame             [cf_start,]  -- object's starting coordinates (default: current coordinates)
                Vector3            [position,]  -- object's new position (default: current position)
                CFrame:Quaternion  [rotation,]  -- object's new rotation (default: current rotation)
                int                [frames,]  -- duration in frames (default: 10)
                LuaFunction        [easing,]  -- easing function (default: linear)
                int                [p1,]      -- parameter for elastic or back type tweening
                int                [p2]       -- parameter for elastic tweening
        ); -- /!\ position and rotation cannot be nil at the same time!

        e.g.    _G.anim.slerp.new(
                        Workspace.Part,
                        nil, -- starts from current coordinates
                        Vector3.new(0, 10, 0),
                        _G.anim.quat_new(0, 270, 0),
                        100,
                        _G.anim.easing.inOutBounce
                );

        ________


8.      _G.anim.slerp.running(
                Instance [object]
        ); -- returns boolean: true whether there's any running slerp on 'object', otherwise false
                               if object == nil, checks for any running thread

        e.g.    _G.anim.slerp.running(Workspace.Part);

        ________


9.      _G.anim.slerp.stop(
                Instance [object]
        ); -- stops every running slerp on 'object' and deletes every pending threads of this object
              if object == nil, stops every running slerp and deletes every pending threads

        e.g.    _G.anim.slerp.stop(Workspace.Part);

        ________


10.     _G.anim.easing.[in, out, inOut, outIn] + [tweener's name](theta, beginning, change, duration)

        e.g.    _G.anim.easing

        To overview all the available easing functions: http://hosted.zeh.com.br/tweener/docs/en-us/misc/transitions.html
        Don't forget to press the 'complete' button on the top right-hand corner to see all the available curves!

________


Known issue:

        Using the slerp function: if you enter a starting cframe relative to 'object', in case the thread
        goes into the pending list, object's old cframe will be retrieved instead of its new position.