# PBRT-Learning
A group of GCLers re-implementing PBRT for code practice.

## Raytracing in one weekend

### Question

- ~~what is **Shadow Acne**~~

- in focus distance (ch12),  what the reason that:

  vertical = focus_dist * viewport_height * v;

  lower_left_corner = origin - horizontal/2 - vertical/2 - focus_dist*w;
  
 - checker_texture::auto sines = sin(10 * p.x()) * sin(10 * p.y()) * sin(10 * p.z());

 - perlin noise

 - **[stb_image](https://github.com/nothings/stb)**  `data = stbi_load() successfully, but value()::data=nullptr.`

 - rotate_y : public hittable:: `rec.set_face_normal(rotated_r, normal); ` is error? and normal-->rec.normal

 - volume: `hit_distance`

### Notes

- **focal distance: **the distance between the projection point and the image plane.
- **focus distance:** the distance between the projection point and the plane where everything is in perfect focus the *focus distance*. 

### Output

- RayTracing in1weekend

  <img src=scene/raytracing_in1weekend/Final-Scene.png width="50%" height="50%">



- RayTracing in the next week

  <img src=scene/RayTracing_next_week/Final_Scene.png width="50%" height="50%">
  
  

- RayTracin TheRestOfYourLife

  <img src=scene/raytracing_the_rest/final_scene.png width="50%" height="50%">



---

## Physically Based Rendering

### Chapter 05

#### **Radiometry**

<img src=trivial/radiometric.png width="70%" height="70%">

1. Basic quantities

   - **Energy:** A photon at wavelength $\lambda$ carries energy
     $$
     Q=\frac{hc}{\lambda}
     $$

   - **Flux:** the limit of differential energy per differential time
     $$
     \Phi=\lim_{\Delta t\rightarrow0} \frac{\Delta Q}{\Delta t}=\frac{dQ}{dt}
     $$

   - **Irradiance and Radiant Exitance:** flux per area

     This quantity is either irradiance(E), **the area density of flux arriving at a surface**, 

     or radiant exitance(M), **the area density of flux leaving a surface.**

     Define: the limit of differential power per differential area at a point p:
     $$
     E(p)=\lim_{\Delta A \rightarrow 0} \frac{\Delta \Phi(p)}{\Delta A}=\frac{d\Phi(p)}{dA}
     $$

     $$
     \int_{A}E(p) d A
     $$

   - **Intensity:**  the limit of a differential cone of directions:
     $$
     I=\lim_{\Delta \omega \rightarrow 0} \frac{\Delta \Phi}{\Delta \omega}=\frac{d\Phi}{d\omega}
        \\
        \Phi=\int_{\Omega} I(\omega)d\omega
     $$
     
   
    Intensity describes the directional distribution of light, but it is only meaningful for point light sources.
   
   - **Radiance:** differential power per differential area **per differential solid angle.**
   
   - $$
     L(p,\omega)=\lim_{\Delta \omega \rightarrow 0} \frac{\Delta E_\omega(p)}{\Delta \omega}=\frac{dE_\omega(p)}{d\omega}
     $$
   
     $E_\omega$: irradiance at the surface that is **perpendicular to the direction** $\omega$
   
     Radiance is the flux density **per unit area, per unit solid angle.**
     $$
     L=\frac{d\Phi}{d\omega dA^{\bot}}
     $$
     <img src=trivial/Radiance.png width="60%" height="50%">

#### Render

$$
E(p,\pmb n)=\int_{\Omega} L_{i}(p,\omega)|cos\theta|d\omega
$$

- **Integrals over spherical coordinates**
  $$
  d\omega=sin\theta \ d\theta \ d\phi\\
  \\
  \begin{align}
  E(p,\pmb n)&=\int_{H^2(\pmb n)} L_{i}(p,\omega)|cos\theta|d\omega\\
  &=\int_0^{2 \pi} \int_0^{\pi/2} L_i(p,\theta,\phi) \ cos\theta\ sin\theta \ d\theta \ d\phi
  \end{align}
  $$
  
- **Integrals over area**
  $$
  d\omega=\frac{dA\ cos\theta}{r^2}\\
  \\
  E(p,\pmb n)=\int_A L\ cos\theta_i\frac{cos\theta_o \ dA }{r^2}
  $$
  <img src=trivial/integral_area.png width="60%" height="50%">

#### Surface reflection

- **BRDF**

  <img src=trivial/brdf.png width="60%" height="50%">

  - the differential irradiance at p is

  $$
  dE(p,\omega_i)=L_i(p,\omega_i)\ cos\theta_i \ d\omega_i
  $$

  - Because of the linearity assumption from geometric optics, the reflecteddifferential radiance is proportional to the irradiance

  $$
  dL_o(p,\omega_o) \propto dE(p,\omega_i)
  $$

  - The constant of proportionality defines the surface’s BRDF for the particular pair of directions $\omega_i$ and $\omega_o$:
    $$
    f_r(p,w_o,\omega_i)=\frac{dL_o(p,\omega_o)}{dE(p,\omega_i)}=\frac{dL_o(p,\omega_o)}{L_i(p,\omega_i)\ cos\theta_i \ d\omega_i}
    $$

- Equation
  $$
  dL_o(p,\omega_o)=f_r(p,w_o,\omega_i)L_i(p,\omega_i)\ |cos\theta_i| \ d\omega_i \\
  \\
  L_o(p,\omega_o)=\int_{\delta^2}f_r(p,w_o,\omega_i)L_i(p,\omega_i)\ |cos\theta_i| \ d\omega_i
  $$

- BSSDF

  <img src=trivial/bssdf_pic.png width="60%" height="50%">

<img src=trivial/bssrdf.png width="80%" height="50%">



### Chapter 02 Geometry and transformations

#### vectors

1. **basic methods:**

   +, -, *, dot, cross, length, normalize, permute, ...

2. Declarations:

   ```c++
   template<typename T> class Vector3{
       public:
       	T operator[](int i) const;
       	T &operator[](int i);
       	Vector3<T> operator+(const Vector3<T> &v) const {
   			return Vector3(x + v.x, y + v.y, z + v.z);
   		}
   		Vector3<T>& operator+=(const Vector3<T> &v) {
   			x += v.x; y += v.y; z += v.z;
   			return *this;
   		}
       ...
           
       public:
       	T x,y,z;
   };
   ```

   ```c++
   typedef Vector2<Float> Vector2f;
   typedef Vector2<int> Vector2i;
   typedef Vector3<Float> Vector3f;
   typedef Vector3<int> Vector3i;
   ```

   

#### points

1. **basic methods:**

   +, -(point), -(vector), distance, Lerp, Permute, operator Vector3<U>(), 



#### Normals

1. **basic methods:**

   Normal and vector are similar, but a normal cannot be **added to a point**,  and one cannot take the **cross product** of two normals

   Faceforward(Normal, Vector); 



#### Rays

1. Form:

$$
r(t)=o+t\pmb d
$$

2. Declarations:

   ```c++
   class Ray {
   	public:
       	Point3f operator()(Float t) const { return o + d * t ; }
       //methods
       
       public:
           Point3f o;
           Vector3f d;
       
       	mutable Float tMax;//it can be changed even if the Ray that contains it is const
       	Float time;
       	const Medium *medium;//each ray records the medium containing its origin. 
   	
   };
   ```



`what is the function of RayDifferential?`

#### Bounding boxes

1. axis-aligned

   - **representation:** the two point **pMin, pMax.**
   - **Methods:** Union, Intersect, bool Overlaps, bool Inside(point, bound), ...

2. Declarations:

   ```c++
   template <typename T> class Bounds3 {
       public:
       	Bounds3(const Point3<T> &p) : pMin(p), pMax(p) { }											Bounds3(const Point3<T> &p1, const Point3<T> &p2);										
       public:
      		Point3<T> pMin, pMax;
       
   };
   ```

   



#### Transformations

1. Homogeneous coordinates
   $$
   (x,y,z,\omega)=(\frac x\omega,\frac y\omega,\frac z\omega)
   $$
   Transformation matrix $M_{4*4}$

2. Basic operations

   - Translations
     $$
     T(\triangle x,\triangle y,\triangle z)=
     \left(
     \begin{array}{l}
     
     1 &0&0&\triangle x\\
     0&1&0&\triangle y\\
     0&0&1&\triangle z\\
     0&0&0&1\\
     \end{array}
     \right)
     $$

   - Scaling
     $$
     S( x,y,z)=\left(\begin{array}{l}
     x &0&0&0\\
     0&y&0&0\\
     0&0&z&0\\
     0&0&0&1\\
     \end{array}\right)
     $$

   - Axis Rotations
     $$
     \pmb R_x(\theta)=\left(\begin{array}{l}
     1 &0&0&0\\
     0&\cos \theta&-\sin \theta&0\\
     0&\sin \theta&\cos \theta&0\\
     0&0&0&1\\
     \end{array}\right)
     $$
     and $\pmb R_y, \pmb R_z$

   - Rotation around an arbitrary axis

     1. **map the given axis to the z-axis.**

     2. rotation in the new system.

     3. recover the origin coords.

        <img src=trivial/1.png width="60%" height="60%">

   - The look-at transformation
   
     The look-at construction  gives a transformation between camera space and world space 
     
     cameraToWorld: **map coordinate system to another.**
     
     <img src=trivial/9.png width="60%" height="60%">
   
3. Declarations:

   ```c++
   class Transform{
       public:
       	Transform Translate(const Vector3f &delta) {
               Matrix4x4 m(1, 0, 0, delta.x,
                           0, 1, 0, delta.y,
                           0, 0, 1, delta.z,
                           0, 0, 0, 1);
               Matrix4x4 minv(1, 0, 0, -delta.x,
                           0, 1, 0, -delta.y,
                           0, 0, 1, -delta.z,
                           0, 0, 0, 1);
               return Transform(m, minv);
           }	
       
   		Transform Scale(Float x, Float y, Float z);
           Transform Rotate(Float theta, const Vector3f &axis); 
       	Transform Rotate(Float theta, const Vector3f &axis);
       	Transform LookAt(const Point3f &pos, const Point3f &look,
                            const Vector3f &up) ;//camera to world
       
       private:
      		Matrix4x4 m, mInv;
       
   };
   ```

   



#### Applying transformations

1. Points and Vectors

    straightforward: $Mp, Mv$

2. Normals

   require special treatment:

   $map(n)=n^*=Sn, S=(M^{-1})^T$

3. Rays

    transforming the constituent **origin and direction**.

4. Bounding boxes

    transform all eight of its corner vertices and then compute a new bounding box.

5. Declarations:

    ```c++
    //Point3f p = ...;
    //Transform T = ...;
    //Point3f pNew = T(p);
    
    //points
    template <typename T> inline Point3<T>
    Transform::operator()(const Point3<T> &p) const {
        T x = p.x, y = p.y, z = p.z;
        T xp = m.m[0][0]*x + m.m[0][1]*y + m.m[0][2]*z + m.m[0][3];
        T yp = m.m[1][0]*x + m.m[1][1]*y + m.m[1][2]*z + m.m[1][3];
        T zp = m.m[2][0]*x + m.m[2][1]*y + m.m[2][2]*z + m.m[2][3];
        T wp = m.m[3][0]*x + m.m[3][1]*y + m.m[3][2]*z + m.m[3][3];
        if (wp == 1) return Point3<T>(xp, yp, zp);
        else return Point3<T>(xp, yp, zp) / wp;
    }
    
    //vectors
    template <typename T> inline Vector3<T>
    Transform::operator()(const Vector3<T> &v) const;
    
    //normals
    template <typename T> inline Normal3<T>
    Transform::operator()(const Normal3<T> &n) const {
        T x = n.x, y = n.y, z = n.z;
        return Normal3<T>(mInv.m[0][0]*x + mInv.m[1][0]*y + mInv.m[2][0]*z,
        mInv.m[0][1]*x + mInv.m[1][1]*y + mInv.m[2][1]*z,
        mInv.m[0][2]*x + mInv.m[1][2]*y + mInv.m[2][2]*z);
    }
    
    //rays
    inline Ray Transform::operator()(const Ray &r) const {
        Vector3f oError;
        Point3f o = (*this)(r.o, &oError);
        Vector3f d = (*this)(r.d);
        //<Offset ray origin to edge of error bounds and compute tMax>
        return Ray(o, d, tMax, r.time, r.medium);
    }
    
    //bounding boxes
    Bounds3f Transform::operator()(const Bounds3f &b) const;
    ```

    

    Common methods:

    ```c++
    Transform Transform::operator*(const Transform &t2) const {
        return Transform(Matrix4x4::Mul(m, t2.m),
        Matrix4x4::Mul(t2.mInv, mInv));
    }
    
    bool Transform::SwapsHandedness() const {
        Float det =
        m.m[0][0] * (m.m[1][1] * m.m[2][2] - m.m[1][2] * m.m[2][1]) -
        m.m[0][1] * (m.m[1][0] * m.m[2][2] - m.m[1][2] * m.m[2][0]) +
        m.m[0][2] * (m.m[1][0] * m.m[2][1] - m.m[1][1] * m.m[2][0]);
        return det < 0;
    }
    ```

    



#### Animating transformations

**interpolate transformations** defined by keyframe matrices based on *matrix decomposition*
$$
M=SRT
$$

- Interpolation of translation and scale: **linear interpolation** 
- Interpolation of  rotations: based on **quaternions.**

**matrix decomposition**

1. T: $T(0:2,3)=M(0:2,3)$
2. R : $M_{i+1}=\frac 12(M_i+(M_i^t)^{-1})$
3. S: $S=R^{-1}M$



#### Interactions

1. Interaction

   - ray intersection: 

     p, normal,scattering media,...

   - Declarations:

     ```c++
     struct Interaction {
     	public:
             Interaction(const Point3f &p, const Normal3f &n, const Vector3f &pError,
         const Vector3f &wo, Float time,
         const MediumInterface &mediumInterface)
         : p(p), time(time), pError(pError), wo(wo), n(n),
         mediumInterface(mediumInterface) { }
         
         public:
         Point3f p;
     	Float time;
         Vector3f pError;// managing floating-point error 
         Vector3f wo;// the negative ray direction
         Normal3f n;
         MediumInterface mediumInterface;
     };
     ```

     

2. surface interaction

   - Variables:
   
     uv, dpdu, dpdv, dndu, dndv, shading,...
   
   - Declarations:
   
     ```c++
     class SurfaceInteraction : public Interaction {
         public:
             SurfaceInteraction(const Point3f &p,
                 const Vector3f &pError, const Point2f &uv, const Vector3f &wo,
                 const Vector3f &dpdu, const Vector3f &dpdv,
                 const Normal3f &dndu, const Normal3f &dndv,
                 Float time, const Shape *shape)
             : Interaction(p, Normal3f(Normalize(Cross(dpdu, dpdv))), pError, wo,
                         time, nullptr),
                 uv(uv), dpdu(dpdu), dpdv(dpdv), dndu(dndu), dndv(dndv),
                 shape(shape) {
                     shading.n = n;
                     shading.dpdu = dpdu;
                     shading.dpdv = dpdv;
                     shading.dndu = dndu;
                     shading.dndv = dndv;
                     //Adjust normal based on orientation and handedness119
                     if (shape && (shape->reverseOrientation ^
                         shape->transformSwapsHandedness)) {
                             n * = - 1 ;
                             shading.n *= -1;
                         }
                 }
         
         	void SetShadingGeometry(const Vector3f &dpdus,
                             const Vector3f &dpdvs, const Normal3f &dndus,
                             const Normal3f &dndvs, bool orientationIsAuthoritative) {
                             //Compute shading.n for SurfaceInteraction
                 				shading.n = Normalize((Normal3f)Cross(dpdus, dpdvs));
                                 if (shape && (shape->reverseOrientation ^
                                             shape->transformSwapsHandedness))
                                         shading.n = -shading.n;
                                 if (orientationIsAuthoritative)
                                     n = Faceforward(n, shading.n);
                                 else
                                 shading.n = Faceforward(shading.n, n);
                             //Initialize shading partial derivative values
                 				shading.dpdu = dpdus;
                                 shading.dpdv = dpdvs;
                                 shading.dndu = dndus;
                                 shading.dndv = dndvs;
                             }
         
         public:
             Point2f uv;
             Vector3f dpdu, dpdv;
             Normal3f dndu, dndv;
             const Shape *shape = nullptr;
         	//shading
             struct {
                 Normal3f n;
                 Vector3f dpdu, dpdv;
                 Normal3f dndu, dndv;
             } shading;
     };
     
     ```
     
     ```c++
     //Transform Method Definitions + ≡
     SurfaceInteraction
     Transform::operator()(const SurfaceInteraction &si) const {
         SurfaceInteraction ret;
         //Transform p and pError in SurfaceInteraction
         //Transform remaining members of SurfaceInteraction
         return ret;
     }
     ```
     
     



---



### Chapter 03 Shapes

#### Basic interface

- ObjectBound()
- WorldBound()
- Intersect(const Ray &ray, Float *tHit,
  SurfaceInteraction *isect, bool testAlphaTexture = true)
- bool IntersectP(const Ray &ray, testAlphaTexture = true)
- Area()

Declarations:

```c++
class Shape{
    public:
    	Shape::Shape(const Transform *ObjectToWorld,
            const Transform *WorldToObject, bool reverseOrientation)
            : ObjectToWorld(ObjectToWorld), WorldToObject(WorldToObject),
            reverseOrientation(reverseOrientation),
            transformSwapsHandedness(ObjectToWorld->SwapsHandedness()) {}
    
    Bounds3f Shape::WorldBound() const {
        return (*ObjectToWorld)(ObjectBound());
    }
    //interface
    public:
    	virtual Bounds3f ObjectBound() const = 0;
        virtual bool Intersect(const Ray &ray, Float *tHit,
                    SurfaceInteraction *isect, bool testAlphaTexture = true) const = 0;
    	//a predicate function
    	virtual bool IntersectP(const Ray &ray,
        bool testAlphaTexture = true) const {
            Float tHit = ray.tMax;
            SurfaceInteraction isect;
            return Intersect(ray, &tHit, &isect, testAlphaTexture);
        }
        virtual Float Area() const = 0;
    public:
    	const Transform *ObjectToWorld, *WorldToObject;
        const bool reverseOrientation;
        const bool transformSwapsHandedness;
}
```

Bounds IntersectP:

```c++
template <typename T>
inline bool Bounds3<T>::IntersectP(const Ray &ray, Float *hitt0,
    Float *hitt1) const {
    Float t0 = 0, t1 = ray.tMax;
    for ( int i = 0 ; i < 3 ; ++i) {
    	Float invRayDir = 1 / ray.d[i];
        Float tNear = (pMin[i] - ray.o[i]) * invRayDir;
        Float tFar = (pMax[i] - ray.o[i]) * invRayDir;
        
        if (tNear > tFar) std::swap(tNear, tFar);
        //UpdatetFarto ensure robust ray–bounds intersection
        t0 = tNear > t0 ? tNear : t0;
        t1 = tFar < t1 ? tFar : t1;
        if (t0 > t1) return false;
    }
    if (hitt0) *hitt0 = t0;
    if (hitt1) *hitt1 = t1;
    return true;
}

//accelerate
template <typename T>
inline bool Bounds3<T>::IntersectP(const Ray &ray, const Vector3f &invDir,
const int dirIsNeg[3]) const {...}
```



#### Spheres

<img src=trivial/10.png width="40%" height="40%">

- basement

$$
\begin{array}{l}
x=r\sin\theta\cos\phi \\
y=r\sin\theta\sin\phi \\
z=r\cos\theta \\
\phi=u\phi_{max}\\
\theta=\theta_{min}+v(\theta_{max}-\theta_{min})
\end{array}
$$

$$
\begin{align}
&\frac{\partial p}{\partial u}=(-\phi_{max}y,\phi_{max}x,0)\\
&\frac{\partial p}{\partial v}=(\theta_{max}-\theta_{min})(z\cos\phi,z\sin\phi,-r\sin\theta)\\
\end{align}
$$

- Weingarten equations

  <img src=trivial/3.png width="60%" height="60%">

- Area

  $$
  A=\phi_{max}r(z_{max}-z_{min})
  $$
  



#### Cylinders

#### Disks

#### Other quadrics

1. Cones
2. Paraboloids
3. Hyperboloids



#### Triangle meshes

1. Triangle intersection

   Methods: Transform the ray such that **its origin is at (0 0 0)** and such that **its direction is along the +z axis.**

   `why not intersect directly?`

   - dpdu, dpdv

     <img src=trivial/4.png width="100%" height="100%">

   - pHit , u,v 

     using  the ***barycentric coordinates***

     <img src=trivial/5.png width="100%" height="100%">

2. Declarations:

   - TriangleMesh:

     ```c++
        struct TriangleMesh {
        	public:
                TriangleMesh(const Transform &ObjectToWorld,
                        int nTriangles, const int *vertexIndices, int nVertices,
                        const Point3f *P, const Vector3f *S, const Normal3f *N,
                        const Point2f *UV,
                        const std::shared_ptr<Texture<Float>> &alphaMask)
                        : nTriangles(nTriangles), nVertices(nVertices),
                        vertexIndices(vertexIndices, vertexIndices + 3 * nTriangles),
                        alphaMask(alphaMask) {
                            //Transform mesh vertices to world space
                            p.reset(new Point3f[nVertices]);
                            for ( int i = 0 ; i < nVertices; ++i)
                                p[i] = ObjectToWorld(P[i]);
                            //Copy UV,N, and S vertex data, if present
                        }
            
            		
            public:
            	const int nTriangles, nVertices;
                std::vector<int> vertexIndices;
                std::unique_ptr<Point3f[]> p;
                std::unique_ptr<Normal3f[]> n;
                std::unique_ptr<Vector3f[]> s;//tangents
                std::unique_ptr<Point2f[]> uv;
                std::shared_ptr<Texture<Float>> alphaMask;
        };
     ```
   
   - Triangle：
   
     ```c++
     class Triangle : public Shape {
         public:
         	Triangle(const Transform *ObjectToWorld, const Transform *WorldToObject,
                     bool reverseOrientation,
                     const std::shared_ptr<TriangleMesh> &mesh, int triNumber)
                     : Shape(ObjectToWorld, WorldToObject, reverseOrientation),
                     mesh(mesh) {
                     v = &mesh->vertexIndices[3 * triNumber];//the triNumber triagnle in mesh
                     }
         
         	//create triangle mesh
         	std::vector<std::shared_ptr<Shape>> CreateTriangleMesh(
                         const Transform *ObjectToWorld, const Transform *WorldToObject,
                         bool reverseOrientation, int nTriangles,
                         const int *vertexIndices, int nVertices, const Point3f *p,
                         const Vector3f *s, const Normal3f *n, const Point2f *uv,
                         const std::shared_ptr<Texture<Float>> &alphaMask) {
                      std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>(
                         *ObjectToWorld, nTriangles, vertexIndices, nVertices, p, s, n, uv,
                         alphaMask);
                         std::vector<std::shared_ptr<Shape>> tris;
                         for ( int i = 0 ; i < nTriangles; ++i)
                             tris.push_back(std::make_shared<Triangle>(ObjectToWorld,
                             WorldToObject, reverseOrientation, mesh, i));
                         return tris;
                         }
         
         	Bounds3f ObjectBound() const {
                     const Point3f &p0 = mesh->p[v[0]];
                     const Point3f &p1 = mesh->p[v[1]];
                     const Point3f &p2 = mesh->p[v[2]];
                     return Union(Bounds3f((*WorldToObject)(p0), (*WorldToObject)(p1)),
                     (*WorldToObject)(p2));
                     }
         	Bounds3f WorldBound() const {
                 const Point3f &p0 = mesh->p[v[0]];
                 const Point3f &p1 = mesh->p[v[1]];
                 const Point3f &p2 = mesh->p[v[2]];
                 return Union(Bounds3f(p0, p1), p2);
                 }
         
         	bool Triangle::Intersect(const Ray &ray, Float *tHit,
                 SurfaceInteraction *isect, bool testAlphaTexture) const {
                		...
                 return true;
                 }
         private:
         	std::shared_ptr<TriangleMesh> mesh;
             const int *v;//triangle_indices:v[0],v[1],v[2]
     };
     ```
     
     
     
     

#### Curves

be used to represent thin shapes for modeling fine geometry like hair, fur, or fields of grass.

#### Subdivision surfaces

#### Managing rounding error





---

### Chapter 04 Primitives and intersection acceleration

#### Primitive interface

1. Declarations:

   ```c++
   class Primitive {
       public:
       	virtual Bounds3f WorldBound() const = 0;
       	//updating Ray::tMax with this value if an intersection is found.
       	virtual bool Intersect(const Ray &r, SurfaceInteraction *) const = 0;
           virtual bool IntersectP(const Ray &r) const = 0;
   		//If the primitive is not emissive, return nullptr
           virtual const AreaLight *GetAreaLight() const = 0;
       	//returns a pointer to the material instance assigned to the primitive.
           virtual const Material *GetMaterial() const = 0;
       	
       	virtual void ComputeScatteringFunctions(SurfaceInteraction *isect,
                       MemoryArena &arena, TransportMode mode,
                       bool allowMultipleLobes) const = 0;
   };
   ```

2. <SurfaceInteraction Public Data> + ≡

   ```c++
   const Primitive *primitive = nullptr;
   //passed to ComputeScatteringFunctions().
   BSDF *bsdf = nullptr;
   BSSRDF *bssrdf = nullptr;
   ```



#### Geometric primitives

The `GeometricPrimitive` class represents a single shape (e.g., a sphere) in the scene.

1. Declarations:

   ```c++
   class GeometricPrimitive : public Primitive {
       public:
       	bool Intersect(const Ray &r,
                   SurfaceInteraction *isect) const {
                   Float tHit;
                   if (!shape->Intersect(r, &tHit, isect))
                       return false;
                   r.tMax = tHit;
                   isect->primitive = this;
                   //Initialize SurfaceInteraction::mediumInterface after Shape intersection
                   return true;
                   }
       
           void GeometricPrimitive::ComputeScatteringFunctions(
                   SurfaceInteraction *isect, MemoryArena &arena, TransportMode mode,
                   bool allowMultipleLobes) const {
                   if (material)
                       material->ComputeScatteringFunctions(isect, arena, mode,
                       allowMultipleLobes);
                   }
       
       private:
           std::shared_ptr<Shape> shape;
           std::shared_ptr<Material> material;
           //this pointer is set to nullptr if the primitive does not emit light
           std::shared_ptr<AreaLight> areaLight;
           MediumInterface mediumInterface;
   };
   ```



#### TransformedPrimitive

Object instancing and Animated primitives.

`TransformedPrimitive` holds a single `Primitive` and also includes an `AnimatedTransform` that is injected in between the underlying primitive and its representation in the scene.

1. Declarations:

   ```c++
   class TransformedPrimitive : public Primitive {
       public:
           TransformedPrimitive(std::shared_ptr<Primitive> &primitive,
               const AnimatedTransform &PrimitiveToWorld)
               : primitive(primitive), PrimitiveToWorld(PrimitiveToWorld) { }
       
       	bool TransformedPrimitive::Intersect(const Ray &r,
                   SurfaceInteraction *isect) const {
                   Transform InterpolatedPrimToWorld;
                   PrimitiveToWorld.Interpolate(r.time, &InterpolatedPrimToWorld);
                   Ray ray = Inverse(InterpolatedPrimToWorld)(r);
                   if (!primitive->Intersect(ray, isect))
                       return false;
                   r.tMax = ray.tMax;
                   //Transform instance’s intersection data to world space
               	if (!InterpolatedPrimToWorld.IsIdentity())
                   *isect = InterpolatedPrimToWorld(*isect);
                   return true;
                   }
       	Bounds3f WorldBound() const {
               return PrimitiveToWorld.MotionBounds(primitive->WorldBound());
           }
       private:
       	std::shared_ptr<Primitive> primitive;
           const AnimatedTransform PrimitiveToWorld;
   };
   ```



#### Aggregates

The `Aggregate` class provides an interface for **grouping multiple Primitive** **objects** together. 

1. Declarations:

   ```c++
   class Aggregate : public Primitive {
       public:
       ...
   };
   ```

   

#### Bounding volume hierarchies

 approaches for ray intersection acceleration based on **primitive subdivision**.

1. Partitioned Algorithms

   - Simple Methods:
     - Split_Middle
     - Split_Equal_Counts
     
   - Other Algorithms
     - Surface area heuristic(SAH)
     
       > 1. rays will incur the cost in splitting the region:
       >
       > $$
       > c(A,B)=t_{trav}+P_A\sum_{i=1}^{N_A}t_{isect}(a_i)+P_B\sum_{i=1}^{N_B}t_{isect}(b_i)\\
       > P(A|C)=\frac{S_A}{S_C}
       > $$
       >
       > 2.  divide the range along the axis into a small number of buckets of equal extent, and then  find the partition with minimum cost.
       >
       >    <img src=trivial/6.png width="60%" height="60%">
       >
       > `why t_trac=1/8 and t_isect=1 `
     
     - Linear bounding volume hierarchies

2. Steps:

   - First, bounding information about each primitive is computed and stored in an array.

   - Next, the tree is built using the algorithm choice.

     > We select the axis associated with the **largest extent when projecting the bounding box**
     > **centroid** for the current set of primitives. 

   - Finally, this tree is converted to a more compact (and thus more efficient) **pointerless representation** for
     use during rendering.
     
     >  The final BVH is stored in a linear array in memory
     >
     > <img src=trivial/7.png width="60%" height="60%">

2. Declarations:

   ```c++
   enum class SplitMethod { SAH, HLBVH, Middle, EqualCounts };
   
   struct BVHPrimitiveInfo {
       BVHPrimitiveInfo(size_t primitiveNumber, const Bounds3f &bounds)
       	: primitiveNumber(primitiveNumber), bounds(bounds),
       	centroid(.5f * bounds.pMin + .5f * bounds.pMax) { }
       size_t primitiveNumber;
       Bounds3f bounds;
       Point3f centroid;
   };
   
   struct BVHBuildNode {
       //BVHBuildNode Public Methods
       void InitLeaf(int first, int n, const Bounds3f &b) {
           firstPrimOffset = first;
           nPrimitives = n;
           bounds = b;
           children[0] = children[1] = nullptr;
           }
       
       void InitInterior(int axis, BVHBuildNode *c0, BVHBuildNode *c1) {
           children[0] = c0;
           children[1] = c1;
           bounds = Union(c0->bounds, c1->bounds);
           splitAxis = axis;
           nPrimitives = 0;
           }
       
       Bounds3f bounds;
       BVHBuildNode *children[2];
       int splitAxis, firstPrimOffset, nPrimitives;
   };
   
   class BVHAccel{
       public:
       	BVHAccel(const std::vector<std::shared_ptr<Primitive>> &p,
           int maxPrimsInNode, SplitMethod splitMethod)
           : maxPrimsInNode(std::min(255, maxPrimsInNode)), primitives(p),
           splitMethod(splitMethod) {
           if (primitives.size() == 0)
               return;
           //Build BVH from primitives
               //first
           std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
           for (size_t i = 0 ; i < primitives.size(); ++i)
               primitiveInfo[i] = { i , primitives[i]->WorldBound() };
               
               //next
               MemoryArena arena(1024 * 1024);
               int totalNodes = 0;
               std::vector<std::shared_ptr<Primitive>> orderedPrims;
               BVHBuildNode *root;
               if (splitMethod == SplitMethod::HLBVH)
                   root = HLBVHBuild(arena, primitiveInfo, &totalNodes, orderedPrims);
               else
                   root = recursiveBuild(arena, primitiveInfo, 0, primitives.size(),
                   &totalNodes, orderedPrims);
               primitives.swap(orderedPrims);
           }
       
       	BVHBuildNode *BVHAccel::recursiveBuild(MemoryArena &arena,
               std::vector<BVHPrimitiveInfo> &primitiveInfo, int start,
               int end, int *totalNodes,
               std::vector<std::shared_ptr<Primitive>> &orderedPrims) {
                   BVHBuildNode *node = arena.Alloc<BVHBuildNode>();
                   (*totalNodes)++;
                   //Compute bounds of all primitives in BVH node
                   int nPrimitives = end - start;
                   if (nPrimitives == 1) {
                   	//Create leaf BVHBuildNode 
                       int firstPrimOffset = orderedPrims.size();
                       for (int i = start; i < end; ++i) {
                           int primNum = primitiveInfo[i].primitiveNumber;
                           orderedPrims.push_back(primitives[primNum]);
                       }
                       node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
                       return node;
                   } else {
                       //Compute bound of primitive centroids, choose split dimension dim
                       Bounds3f centroidBounds;
                       for (int i = start; i < end; ++i)
                           centroidBounds = Union(centroidBounds, primitiveInfo[i].centroid);
                       int dim = centroidBounds.MaximumExtent();
                       //Partition primitives into two sets and build children
                       int mid = (start + end) / 2;
                       if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
                           //Create leaf BVHBuildNode
                       } else {
                       //Partition primitives based on splitMethod
                           node->InitInterior(dim,
                               recursiveBuild(arena, primitiveInfo, start, mid,
                               totalNodes, orderedPrims),
                               recursiveBuild(arena, primitiveInfo, mid, end,
                               totalNodes, orderedPrims));
                       }
                   }
                   return node;
               }
       
       private:
       	const int maxPrimsInNode;
           const SplitMethod splitMethod;
           std::vector<std::shared_ptr<Primitive>> primitives;
   };
   
   ```
   



#### KD-Tree

ray intersection acceleration based on **spatial partition.**

 A kd-tree simply restricts the splitting plane to be perpendicular to one of the coordinate axes; 

<img src=trivial/8.png width="60%" height="60%">



### Reflcetion Model

https://www.zhihu.com/search?type=content&q=%E5%BE%AE%E8%A1%A8%E9%9D%A2%E6%A8%A1%E5%9E%8B

