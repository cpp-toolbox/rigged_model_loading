# rigged_model_loading

As of right now the rigged model loading only seems to support `.fbx` so make sure you export to this format, also when exporting you need to make sure you have specific setting on or else there can be problems. As of right now I only have documentation to do this in blender: 

![image](https://github.com/user-attachments/assets/3b2a745d-437b-46e5-8f8a-90ac8d9cb0f4)

This has to occur because there are various problems with armatures when they get exported to fbx having a scaledown of 100x and also being rotated usually.

## multiple armature support

**Warning: when using multiple armatures you need to make sure that the bone names across all armatures are unique and that the armature names themselves are unique, and follow the naming convention in the next paragraph.**

*Note: In the future we will try and remove these requirements to simplify things.*

You must name your armatures using the following format `X_..._armature` eg) `character_armature` and so on, also animations must follow a similar convention, so for example you have to name your action `character_anim`, this must be done because when exporting an animation with mutiple armatures with the animation it generates a bunch of empty animations for the armatures that were not used during that animation, and we needed a way to extract which animation is for what armature, this also works because when blender exports to fbx the actions get named like this: `character_armature|character_anim`, so that when we see an animation of this form`character_armature|gun_anim` then we don't have to care about it.

## animations with constraints

When creating animations you usally use the child-of constraint to make an object follow another object, such has keeping an object in a characters hand while they move the hand around. Note that managing the transforms of objects that use these constraints can become very combersome quickly when transfering objects between parents, if you work on animations that use a lot of re-parenting I thoroughly recommend the dynamic parent plugin.

Since constraints are a blender specific thing and not specified in regular keyframes which only contain rotation, translation and scale of bones, then when exporting you have to make sure of the following: 

*  bake your actions <b>make sure to keep the cursor at 0 when you bake the animation (I'm not sure why this has to happen yet)</b>
* When baking your actions be sure all of the individual armatures are selected if you don't do this and some of the animations on your armatures start a different moments in time then their animations will not be in sync
* Mka eusre that that you give all armatures some key frames that assimp believes they are animated.

 which pretty much gets rid of the constraints and instead bakes in the translation that would have resulted from that constraint, allowing for the exported format animation to contain all the data required to render the animation, to do this open the bake action dialog and use these settings: 

![image](https://github.com/user-attachments/assets/29d8b148-b5bf-41d7-b3c7-9ff417de1330)


**Note: When using the child of constraint it's important that you set the inverse transform relative to the origin of the scene or else when you parent it might not work correctly**

