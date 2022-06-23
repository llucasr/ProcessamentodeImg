//==============================================================================================
// Originally written in 2016 by Peter Shirley <ptrshrl@gmail.com>
//
// To the extent possible under law, the author(s) have dedicated all copyright and related and
// neighboring rights to this software to the public domain worldwide. This software is
// distributed without any warranty.
//
// You should have received a copy (see file COPYING.txt) of the CC0 Public Domain Dedication
// along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//==============================================================================================

#include "rtweekend.h"

#include "camera.h"
#include "color.h"
#include "hittable_list.h"
#include "material.h"
#include "sphere.h"
#include <fstream>
#include "stb_image.h"
#include "stb_image_write.h"
#include <iostream>


color ray_color(const ray& r, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0,0,0);

    if (world.hit(r, 0.001, infinity, rec)) {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth-1);
        return color(0,0,0);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}


hittable_list random_scene(double px, double py) {

    hittable_list world;


    auto ground_material = make_shared<lambertian>(color(0.3, 0.7, 0.7));
    world.add(make_shared<sphere>(point3(0,-1000,0), 1000, ground_material));

   // for (int a = 0; a < 1; a++) {
    //    for (int b = 0,i=0; b < 1; b++) {


            auto choose_mat = random_double();
            //for(u=0;u<1.0;u=u+0.02){

            point3 center(px,0.2,py);
            std::cout<<px<<")("<<py<<"\n";


          //  point3 center(a + 0.5*random_double(), 0.2, b + 0.9*random_double());

           // center - point3(i+=0.2, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;
                auto albedo = color::random() * color::random();
                sphere_material = make_shared<lambertian>(albedo);
                world.add(make_shared<sphere>(center, 0.2, sphere_material));
/*}
                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                } else {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }*/
       // }
   // }

    /*
    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));
*/
    return world;
}


int main() {

    // Image

    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 1200;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 10;
    const int max_depth = 50;
    int k = 0;
    char nomeimg[11];
    int i = 0;

    double px[202], py[202], u= 0;
    double x[6],y[6];

    x[0]=y[0]=5;
    x[1]=1;y[1]=4;
    x[2]=2;y[2]=4;
    x[3]=1;y[3]=0;



    for(u=0;u<1.0;u=u+0.005,k++){

    px[k] = pow(1-u,3)*x[0]+3*u*pow(1-u,2)*x[1]+3*pow(u,2)*(1-u)*x[2]+pow(u,3)*x[3]    ;
    py[k] = pow(1-u,3)*y[0]+3*u*pow(1-u,2)*y[1]+3*pow(u,2)*(1-u)*y[2]+pow(u,3)*y[3];

    }
    std::cout<<k;
    for(i=0;i<1;i++){

// World

    auto world = random_scene(px[i],py[i]);

    // Camera

    point3 lookfrom(-3,2,1);
    point3 lookat(0,0,-1);
    vec3 vup(0,1,0);
    auto dist_to_focus = 11.0;
    auto aperture = 0.0;


    camera cam(lookfrom, lookat, vup, 90, aspect_ratio, aperture, dist_to_focus);

    // Render

    sprintf(nomeimg,"img%d.ppm",i);

    std::ofstream out(nomeimg);


 out << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    for (int j = image_height-1; j >= 0; --j) {
        std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0,0,0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width-1);
                auto v = (j + random_double()) / (image_height-1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            write_color(out, pixel_color, samples_per_pixel);

        }}
       std::cout<<"\n";

    }


}

