# OctreeImplementation
Octree Implementation in huge 3D static space, to decrese the time bewteen frames in a small engine in comparison with a object collision detection done by check every player with every object.

My solution for this test about a huge 3D static space, taking in count that the most important optimization that must be done is decrease the time between frames, is an Octree implementation.

Results:
I lose time in the initialization of the world and memory use, in exchange of a huge increase in performance of frame rate.

CPU: i7-6700HQ @ 2.60GHz
Platform: x64

                              Initialization Time 	Average time between frames
    Vanilla				3500.86ms               1370.25ms
    With Octree implementation	13712.6ms		0.0386ms
