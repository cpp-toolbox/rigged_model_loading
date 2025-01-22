# rigged_model_loading

As of right now the rigged model loading only seems to support `.fbx` so make sure you export to this format, also when exporting you need to make sure you have specific setting on or else there can be problems. As of right now I only have documentation to do this in blender: 

![image](https://github.com/user-attachments/assets/3b2a745d-437b-46e5-8f8a-90ac8d9cb0f4)

This has to occur because there are various problems with armatures when they get exported to fbx having a scaledown of 100x and also being rotated usually.

## multiple armature support

**Warning: when using multiple armatures you need to make sure that the bone names across all armatures are unique and that the armature names themselves are unique, and follow the naming convention in the next paragraph.**

*Note: In the future we will try and remove these requirements to simplify things.*

You must name your armatures using the following format `X_..._armature` eg) `character_armature` and so on, also animations must follow a similar convention, so for example you have to name your action `character_anim`, this must be done because when exporting an animation with mutiple armatures with the animation it generates a bunch of empty animations for the armatures that were not used during that animation, and we needed a way to extract which animation is for what armature, this also works because when blender exports to fbx the actions get named like this: `character_armature|character_anim`, so that when we see an animation of this form`character_armature|gun_anim` then we don't have to care about it.

## Baking Animations

### Why Bake?

When making animations in blender you have access to tools that will make animating easier, usually these tools are not portable and are specific to blender, thus when you export your animation, those tools will not be able to be used, thus you must "bake" in the tools effect into the file on export or else other programs which try to open the file will not see these effects. A simple example of this would be if you only work in `.blend` files and are using a mirror modifier and you export to `.fbx` or `.obj` clearly that reflect data will not be there and you'd only get half the model, this is why baking in general is important.

### So what data is exported?

From what I understand when you export an animation you get keyframes which only contain rotation, translation and scale of bones

### Interpolation Baking

Blender is the tool for creating animations, the game engine does not, thus when we do fancy interpolation, then that information would not be exported into the fbx file usually, game engines usually just store the keyframes and its up to the program to do the interpolation during the animation. Due to this constraint, it's best that you "bake" this data into the animation while exporting, so that the other program can just linearly interpolate in small steps, which will still be a good approximation of your fancy interpolation in blender.

### Constraint Baking

When creating animations you usally use the child-of constraint to make an object follow another object, such has keeping an object in a characters hand while they move the hand around. This is usually blender specific and these constraints will not be in the exported file, thus we have to bake this too.

**Note**: Managing the transforms of objects that use these constraints can become very combersome quickly when transfering objects between parents, if you work on animations that use a lot of re-parenting I thoroughly recommend the dynamic parent plugin.

### How to Bake

* Keep the timeline cursor at 0 when you bake the animation (I'm not sure why this has to happen yet)</b>
* in object mode, select all your armatures
* go to pose mode
* open the bake action dialog and use these settings: 

![image](https://github.com/user-attachments/assets/29d8b148-b5bf-41d7-b3c7-9ff417de1330)

**Note**:  Make eusre that that you give all armatures some key frames that assimp believes they are animated.


**Note**: When using the child of constraint it's important that you set the inverse transform relative to the origin of the scene or else when you parent it might not work correctly

