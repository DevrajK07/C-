#include "Header.h"


void xacceleration(body *bodies, int current_body, int n, double G)
{
    bodies[current_body].ax = 0;
    for(int i = 0; i < n; i++)
    {
        if(i != current_body)
        {
            double d = sqrt(pow((bodies[i].x - bodies[current_body].x),2) + pow((bodies[i].y - bodies[current_body].y),2));
            bodies[current_body].ax += (G * bodies[i].m / pow(d,3)) * (bodies[i].x - bodies[current_body].x);
        }
    }
}

void yacceleration(body *bodies, int current_body, int n, double G)
{
    bodies[current_body].ay = 0;
    for(int i = 0; i < n; i++)
    {
        if(i != current_body)
        {
            double d = sqrt(pow((bodies[i].x - bodies[current_body].x),2) + pow((bodies[i].y - bodies[current_body].y),2));
            bodies[current_body].ay += (G * bodies[i].m / pow(d,3)) * (bodies[i].y - bodies[current_body].y);
        }
    }
}

//Fourth Order Runge-Kutta Function
void rungeKutta(body *bodies, int n, double t, double dt, double G)
{
    // Defining fourth order Runge Kutta terms
    double k1x, k2x, k3x, k4x; //x displacements
    double k1y, k2y, k3y, k4y; //y displacements
    double k1vx, k2vx, k3vx, k4vx; //x velocities
    double k1vy, k2vy, k3vy, k4vy; //y velocities
    
    // For loop to cycle through the number of bodies
    for(int i=0; i<n; i++)
    {
        //Assumes linear change in acceleration between each time step dt (for more accurate outcome decrease dt)
        //Substitutes the new bodies, current body and number of bodies per iteration of RK4 into acceleration function
        xacceleration(bodies, i, n,G);
        yacceleration(bodies, i, n,G);

        k1vx = bodies[i].ax * dt; //x velocity
        k1vy = bodies[i].ay * dt; //y velocity
        k1x = bodies[i].vx * dt; //x displacement
        k1y = bodies[i].vy * dt; //y displacement

        bodies[i].x += 0.5*k1x;
        bodies[i].y += 0.5*k1y;
        bodies[i].vx += 0.5*k1vx;
        bodies[i].vy += 0.5*k1vy;
        xacceleration(bodies, i, n,G);
        yacceleration(bodies, i, n,G);
        
        //Calculate k2's
        k2vx = (bodies[i].ax) * dt; //x velocity
        k2vy = (bodies[i].ay) * dt; //y velocity
        k2x = (bodies[i].vx)*dt; //x displacement
        k2y = (bodies[i].vy)*dt; //y displacement

        bodies[i].x += -0.5*k1x + 0.5*k2x;
        bodies[i].y += -0.5*k1y + 0.5*k2y;
        bodies[i].vx += -0.5*k1vx + 0.5*k2vx;
        bodies[i].vy += -0.5*k1vy + 0.5*k2vy;
        xacceleration(bodies, i, n,G);
        yacceleration(bodies, i, n,G);

        //Calculate k3's
        k3vx = (bodies[i].ax)* dt; //x velocity
        k3vy = (bodies[i].ay)* dt; //y velocity
        k3x = (bodies[i].vx)*dt; //x displacement
        k3y = (bodies[i].vy)*dt; //y displacement

        bodies[i].x += -0.5*k2x + k3x;
        bodies[i].y += -0.5*k2y + k3y;
        bodies[i].vx += -0.5*k2vx + k3vx;
        bodies[i].vy += -0.5*k2vy + k3vy;
        xacceleration(bodies, i, n,G);
        yacceleration(bodies, i, n,G);
    
        //Calculate k4's
        k4vx = (bodies[i].ax) * dt; //x velocity
        k4vy = (bodies[i].ay) * dt; //y velocity
        k4x = (bodies[i].vx)*dt; //x displacement
        k4y = (bodies[i].vy)*dt;   //y displacement
        
        bodies[i].x += -k3x;
        bodies[i].y += -k3y;
        bodies[i].vx += -k3vx;
        bodies[i].vy += -k3vy;

        //Update positions and velocities
        bodies[i].x += (1.0/6.0)*(k1x + 2*k2x + 2*k3x + k4x);
        bodies[i].y += (1.0/6.0)*(k1y + 2*k2y + 2*k3y + k4y);
        bodies[i].vx += (1.0/6.0)*(k1vx + 2*k2vx + 2*k3vx + k4vx);
        bodies[i].vy += (1.0/6.0)*(k1vy + 2*k2vy + 2*k3vy + k4vy);
    }
}