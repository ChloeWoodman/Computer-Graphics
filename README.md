[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-c66648af7eb3fe8bc4f294546bfd86ef473780cde1dea487d3c4ff354943c9ae.svg)](https://classroom.github.com/online_ide?assignment_repo_id=9876457&assignment_repo_type=AssignmentRepo)
# Computer Graphics

## Where to access work
Each week has their own individual branch that store their labs. When all work is complete these will be merged into main.

## Pulling the repo 
`git clone https://github.com/RyanWestwood/Computer-Graphics.git`  
`git submodule update --init --recursive`

## Building the project 
run `build-scripts/build.bat`

## Weekly updates
Week 1:
Using PBRT and cmake to render a scene, at first I was having problems with the pbrt not being a command, as informed by CommandLine, when I attempted to render a scene (kilerroo due to moodle being down, it was a downloaded scene already given from pbrt). I quickly learned that the killeroo needs to be put inside of the release folder (release because it works quicker than debug) so the pbrt is in the file path for commandLine to use. After I did this, the pbrt outputted the scene. I also noticed if you open the og scene in notepad, you can view the features of the scene, like models, textures and lighting - meaning if you play around with these numbers it can change the output of the scene.
![image](https://user-images.githubusercontent.com/113985493/221229259-a31e304f-0317-4f8c-9ebe-2566c7c67d65.png)
Would display output but GitHub does not like pbrt, I even did a copy of it on the file and pushed but it still doesn't work. 


Week 2:
Week 2 consisted of creating C++ code to write images in low and high dynamic range. I was given a TGA file for low dynamic range and a PPM format for high dynamic range format. I cloned from a git repo that had incomplete code, the code would be built in VS and ran to have PPM output rainbow RGB colours andr TGA to output a red dot. The aim was to make TGA match the PPM image in colours. I noticed a lot of the PPM code can be recycled as it uses RGB values, so after putting the right rendering code in the right places, TGA outputted the same as PPM.
![image](https://user-images.githubusercontent.com/113985493/221229411-4eb4375d-4560-41f5-9f2b-5d8b7a9c0c7b.png)


Week 3:
This was the beginning of line drawing, cloned by a repo containing code by extending the C++ code we learned to write images in TGA format via line drawing. Pixel positions (x, y) could be drew with an increase using for loops and assigned rgb values which means a coloured line could be created. Essentially, interpolating a line by drawing pixels with a start and end point. This can then be used to render wireframes of imported models which creates images.
![image](https://user-images.githubusercontent.com/113985493/221229609-07a2f2d7-b203-4cbc-a196-a1ed1440a1fb.png)


Week 4:
This code closely builds off week 3, once again a clonable repo has been supplied, though the code is very similar to what I made last week. This code draws lines and colours, however can also have a routine to draw triangles that was extended upon - triangle drawn can create lines between vertices for definition, which then can be coded to fill the triangle in with the assigned rgb value. The code has an obj parser and model too which allows for when triangles worked and connected, it could create triangles and shade them of a model mesh.
![image](https://user-images.githubusercontent.com/113985493/221230126-865c30f1-6892-46ab-a8ab-6006311bbc61.png)


Week 5:
This week I did not get a lot of work done, I originally started the week 5 lighting challenge lab but quickly changed my mind to work on my Maya scene for milestones, my Maya scene is an aquatic enviroment with a jellyish as the focal point - I intend to have background miscellaneous items as well to display different textures and materials. I have used the Midjourney bot in the BCU discord to render my image by entering /image prompts. Here is the output references for the scene I intend to do:
![image](https://user-images.githubusercontent.com/113985493/221227312-e2186b6b-285e-4090-b9dc-815284aeec62.png)
![image](https://user-images.githubusercontent.com/113985493/221227336-93426876-79bb-415d-b148-25900a51d142.png)
![image](https://user-images.githubusercontent.com/113985493/221227378-bbd2ab3a-64b9-45c9-99ac-b8d4990b0dbe.png)
![image](https://user-images.githubusercontent.com/113985493/221227413-0947f97b-8f67-4fd1-9f2c-89d9212d4e73.png)
(I did this process a few times so I would not be limited to a singular picture and could take inspiration from all of them)

Week 6:
Just working on game scene in Maya for milestone, see Game scene for more details on what the Maya scene ended up looking at.

Week 7: 
DMT futures.

Week 8:
Finished all requirements for a pass criteria for a rasterizer with a working camera!
![image](https://user-images.githubusercontent.com/113985493/226122177-4bb2472c-2d32-454c-93f9-928740478808.png)

Week 9:
This week I moved onto ray tracing and made a red sphere that is written to a PPM file with lighting and shading :).
![image](https://user-images.githubusercontent.com/113985493/228382009-164fc7e8-023f-48a5-817a-14c981e09374.png)


Week 10:
This week I ray traced 3 spheres each with their own different materials, using metal and lambert with assigned colours that can mirror the other spheres if using a reflective material with an applied fuzz.
![image](https://user-images.githubusercontent.com/113985493/228380610-3d1f6ba0-3348-4ad8-9f31-9b474e294189.png)

Week 11(and week 12)
This week I decided to do the labs from both week 11 and week 12 in one so if I had any issues I would have ample time to be able to get support or debug the issue. Now my ray tracer can render in models and sphere scenes, depending on which scene I make as default. They can render in multiple different materials with a perspective camera and motion blur with point lighting. Note the green hue was set by me, it is not an error.
![image](https://user-images.githubusercontent.com/113985493/233843517-03709948-dc51-4f56-af94-e46d6e9af500.png)
![image](https://user-images.githubusercontent.com/113985493/233844335-07647051-3c3d-4a0d-a0f5-7e6e669d086e.png)


# Game scene
![image](https://user-images.githubusercontent.com/113985493/221852767-3fef8ce7-bc80-4235-a271-a58bcf7aa7de.png)
Using Arnold renderer: 
![image](https://user-images.githubusercontent.com/113985493/221861362-012b863d-1933-4538-836e-14e436bd552a.png)



# References
"Seaweed 2" (https://skfb.ly/6SAyy) by lyningsknallis is licensed under Creative Commons Attribution (http://creativecommons.org/licenses/by/4.0/).

"Coral" (https://skfb.ly/6U9xn) by Lydechaser is licensed under Creative Commons Attribution (http://creativecommons.org/licenses/by/4.0/).

"Lowpoly coral" (https://skfb.ly/o9zJx) by assetfactory is licensed under Creative Commons Attribution (http://creativecommons.org/licenses/by/4.0/).

"Jellyfish" (https://skfb.ly/oDXYN) by ZIRODESIGN is licensed under Sketchfab Standard Licensing (https://sketchfab.com/licenses)

"Rudd Fish" (https://www.turbosquid.com/3d-models/free-obj-model-fish-rudd/603242) by viebi is licensed under TurboSquid Standard 3D Model License (https://blog.turbosquid.com/turbosquid-3d-model-license/)

# Resources
Resources used to help with this code's ray tracer was sourced from Ray Tracing in One Weekend book series: [_Ray Tracing in One Weekend_](https://raytracing.github.io/books/RayTracingInOneWeekend.html)
