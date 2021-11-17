# PBRT-Learning
A group of GCLers re-implementing PBRT for code practice.

## Question

- ~~what is **Shadow Acne**~~

- in focus distance (ch12),  what the reason that:

  vertical = focus_dist * viewport_height * v;

  lower_left_corner = origin - horizontal/2 - vertical/2 - focus_dist*w;
  
 - checker_texture::auto sines = sin(10 * p.x()) * sin(10 * p.y()) * sin(10 * p.z());

 - perlin noise

 - **[stb_image](https://github.com/nothings/stb)**  `data = stbi_load() successfully, but value()::data=nullptr.`

 - rotate_y : public hittable:: `rec.set_face_normal(rotated_r, normal); ` is error? and normal-->rec.normal

 - volume: `hit_distance`

## Notes

- **focal distance: **the distance between the projection point and the image plane.
- **focus distance:** the distance between the projection point and the plane where everything is in perfect focus the *focus distance*. 

## Output

- RayTracing in1weekend

  <img src=scene/raytracing_in1weekend/Final-Scene.png width="50%" height="50%">



- RayTracing in the next week

  <img src=scene/RayTracing_next_week/Final_Scene.png width="50%" height="50%">
  
  

- RayTracin TheRestOfYourLife

  <img src=scene/raytracing_the_rest/final_scene.png width="50%" height="50%">
