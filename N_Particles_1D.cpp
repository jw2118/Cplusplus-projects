#include <iostream>
#include <math.h>
#include <stdlib.h> 
#include <time.h>
using namespace std;
//Define global constants, number of particles, NC2 for convenience
const int N=11;
const int NChoose2=(N+2)*(N+1)/2;

//Particle structure to be used to define several particles
struct Particle{
    float position;
    float velocity;
    float mass;
};

//Assigns velocities
void assignPosVelMass(Particle particles[N+2], float pos[N], float vel[N], float mass[N]);
//Min collision time
void collidetime(Particle particles[N+2], float output[3]);
//This function runs the programme forward
void runtime(Particle particles[N+2], float* time);
//Generating random
void randomPos(float pos[N], float width);
void randomVel(float vel[N], int vmax);

int main(){
    srand (time(NULL));
    //Width of the box
    float width = 20;
    //Total time of the universe
    float time=0;
    //Number of particles, fixed constant neither created nor destroyed
    //Last two particles are the walls
    Particle theboys[N+2];
    //Leftwall
    theboys[N].velocity=0;
    theboys[N].position=0;
    //Rightwall
    theboys[N+1].velocity=0;
    theboys[N+1].position=width;
    float initvel[N];
    float initpos[N];
    float masses[N];
    randomPos(initpos, width);
    randomVel(initvel, 10);
    
    for(int i=0; i<N-1;i++){
        masses[i] = 1;
    }

    assignPosVelMass(theboys, initpos, initvel, masses);
    for(int i=0;i<100000;i++){
        runtime(theboys,&time);
    }

    //Now final print statement to put out velocities
    for(int i=0;i<N;i++){
        if(i== N-1){
            //cout << "Initial velocity " << i << " " << initvel[i] << endl;
            cout << theboys[i].velocity << endl;
        }else{
            cout << theboys[i].velocity << ",";
        }
    }
};

void assignPosVelMass(Particle particles[N+2], float pos[N], float vel[N], float mass[N]){
    for(int i=0; i <N; i++){
        particles[i].position = pos[i];
        particles[i].velocity = vel[i];
        particles[i].mass = mass[i];
    }
};

void collidetime(Particle particles[N+2], float output[3]){
    
    //This is half the number but rounded up
    int K = (N+2)/2+(N+2)%2;
    //Dummy variable to start the loop from
    int L =0;

    //Dummy variable to add 1 to each loop to store values correctly
    int Index=0;
    
    //Variable to store the minimumtime
    int mintime=2000;
    //Only negative collision times correspond to a collision occurring
    float relVel[NChoose2];
    float relPos[NChoose2];
    float timeCollide[NChoose2];
    output[0] = mintime;
    //Loop to get all of the particles. We halve to avoid getting repeats and ignore i=j since these will always give zero
    for(int i=0; i<=K; i++){
        L+=1; //This stops repeats of things like 01 and 10 in the following loop.
        for(int j=L; j<N+2; j++){
            if(i!=j){ //This will just give zero, eliminate i=j
                relVel[Index] = particles[i].velocity - particles[j].velocity; 
                relPos[Index] = particles[i].position - particles[j].position;
                timeCollide[Index] = relPos[Index]/relVel[Index];
                //cout << "Collision time for " << i << j << " " << timeCollide[Index] << endl; //Velocity is zero
                //This if statement determines the minimum time and its index
                if(timeCollide[Index]<0 && abs(timeCollide[Index]) < abs(output[0])){
                    //Store minimum time and the indices of the particles colliding
                    output[0] = timeCollide[Index];
                    output[1] = i;
                    output[2] = j;
                }
                Index+=1;
            }  
        }
    
    }
};

//time is the total time elapsed
void runtime(Particle particles[N+2],float* time){
    //Create arrays to store outputs of wall and particle collisions.
    //First index of these is the time and following ones are the indices of the particles.
    float outputData[N+1]; //outputData has time in its first spot, then all the particles in the following spots
    float particlecollisionData[3];
    //Find min collision and between which particles
    collidetime(particles, particlecollisionData);

    float deltaT = particlecollisionData[0];
    int p1Index = particlecollisionData[1];
    int p2Index = particlecollisionData[2];
    *time -= deltaT; //Collision times negative by convention (not really working)
    //First want to handle the particles all changing position, only need to consider first N since walls don't move
    outputData[0] = *time;
    for(int i=0; i<N; i++){
        particles[i].position = particles[i].position-particles[i].velocity*deltaT; 
        outputData[i+1] =particles[i].position;
    }

    int Index =0;
    int K= (N+2)/2+(N+2)%2;
    int L=0;
    float Pos1Pos2[NChoose2];
    if(abs(deltaT) > 0.0001){
        for(int i=0; i<=K; i++){
            L+=1; //This stops repeats of things like 01 and 10 in the following loop.
            for(int j=L; j<N+2; j++){
                if(i!=j){ //Particles don't self-interact
                    if(abs(particles[i].position - particles[j].position) < 0.0001){ // If particles are colliding
                        if(j>=N){//Colliding with wall 
                            particles[i].velocity = -particles[i].velocity;
                        }else if(i>=N){ //Colliding with wall
                            particles[j].velocity = -particles[j].velocity;
                        }else{
                            float v1 = particles[i].velocity;
                            float v2 = particles[j].velocity;

                            float m1 = particles[i].mass;
                            float m2 = particles[j].mass;
                            //ZMF is Zero Momentum Frame velocity 
                            float vZMF = (m1*v1+m2*v2)/(m1+m2);
                            particles[i].velocity = 2*vZMF - v1; 
                            particles[j].velocity = 2*vZMF - v2;
                        }
                    }
                }           
            }
        }

        //Now print all the outputs
        for(int i=0; i<=N;i++){
            if(i==N){
                cout << outputData[i];
            }else{
                cout << outputData[i] << ", ";
            }
        }
        cout << endl;
        
    };
}
//For some reason double switches occur, tracks collision twice and this results in what appears to be no collision even though it knows.
void randomPos(float pos[N], float width){
    int widthInt = (int)width; //COME FIX LATER
    for(int i=0;i<N;i++){
        pos[i] = rand() % widthInt + float(rand()%10000)/10000;
    }
}

void randomVel(float vel[N], int vmax){
    for(int i=0;i<N;i++){
        vel[i] = rand() % vmax - vmax/2 + float(rand()% 10000)/10000;
    }
}