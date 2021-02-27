#include <iostream>
#include <math.h>
using namespace std;
//By default, positive G is attractive and negative G is repulsive
const int G=1;

//Attraction has fixed strength M=1, omitted. 

class vec2D {
    public:
        float x,y;

        vec2D() {
            x=0.0; y=0.0;
        }

        vec2D( float a, float b) {
            x=a; y=b;
        }
        //Friend allows the class to access things that are normally protected
        friend ostream& operator << (ostream &out, const vec2D &v){
            out << v.x << "," << v.y;
            //out << "First: " << v.x << ", Second: " << v.y;
            return out;
        }

        vec2D operator+ (const vec2D &v) const {
            return vec2D(x + v.x, y + v.y);
        }

        vec2D operator- (const vec2D &v) const {
            return vec2D(x + v.x, y + v.y);
        }
        //Only right multiplication by a scalar works
        vec2D operator* (const float &a) const {
            return vec2D(x*a, y*a);
        }

        vec2D operator/ (const float &a) const {
            return vec2D(x/a, y/a);
        }

        float abs2() const {
            return (x*x + y*y); 
        }
        
        float abs() const {
            return sqrt(abs2());
        }

        float firstComp() const{
            return x;
        }

        float secondComp() const{
            return y;
        }



};

//Planet struct to store position, velocity and mass
struct Planet{
    vec2D position;
    vec2D velocity;
    float mass;
};


//Functions: calculating force from position with arbitrary power law and iterator
vec2D Force(float power, Planet planet);
//Total time, t, iterating step,dt
void stepTime(Planet& planet, float power, float t, float dt);

float keplerEnergy(Planet planet);
float angularMomentum(Planet planet);
float Eccentricity(Planet planet);

int main(){
    //vec2D v1(1,2);
    //vec2D v3(-1,4.5);

    //Initialize time and timestep variables
    float t=0;
    float dt=0.01;
    float endTime = 10000;
    float stepNumber = endTime/dt;

    //Define the initial conditions and properties of planet
    Planet Earth;
    Earth.mass = 1;

    //Iterator step to run the programme forward in time
    //large for loop sets the initial velocity for different values
    for(int j=0; j<5; j++){
        float vy = 0.5 * cos(0.25*j);
        float vx = 0.5 * sin(0.25*j);
        Earth.position = vec2D(4,0);
        Earth.velocity = vec2D(vx,vy);
        cout << "Eccentricity: " <<Eccentricity(Earth) << endl; 
        for(int i=0; i < stepNumber; i++){
            //velocity and poisition plot generating loop
            
            if(i%100 ==0){
                //cout << Earth.velocity << endl;
                cout << Earth.position << endl;
            }
            stepTime(Earth, 2, t, dt);
            t= t+dt;
        }
    }
    
    return 0; 
};


vec2D Force(float power, Planet planet){
    vec2D r = planet.position;
    float R = r.abs();
    float k = -planet.mass*G*1/pow(R,power+1);
    return r*k;
};

//The ampersand means you're currently passing the planet by reference i.e its values will be modified
void stepTime(Planet& planet, float power, float t, float dt){
//Start with a basic Euler method
    planet.velocity= planet.velocity + Force(power, planet)*dt;
    planet.position= planet.position + planet.velocity*dt;
};

float keplerEnergy(Planet planet){
    vec2D r = planet.position;
    float R = r.abs();
    float m = planet.mass;
    float vsquare = planet.velocity.abs2();
    float E = -G*m/R + 0.5*m*vsquare;
    return E;
}

float angularMomentum(Planet planet){
    vec2D r = planet.position;
    vec2D v = planet.velocity;
    float m = planet.mass;
    float x = r.firstComp();
    float vy = v.secondComp();
    return m*vy*x;
}

float Eccentricity(Planet planet){
    float E = keplerEnergy(planet);
    float L = angularMomentum(planet);
    float m = planet.mass;
    float k = 1+2*E*L*L/(G*G*m*m);
    float e = sqrt(k);
    return e;
}