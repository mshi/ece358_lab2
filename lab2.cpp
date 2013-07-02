#include <time.h>
#include <iostream>
#include <string>
#include <cstdlib>
#define DEBUG
#include "Debug.h"
#include "random.h"
#include "simulation.h"

using namespace std;

int main(int argc, char* argv[]) {
  // Seed the RNG
  sgenrand(time(NULL));

  unsigned long numConnected;  // number of computers connected to LAN
  double arrivalRate;  // arrival rate A as packets/second
  double persistence;
  const double lanSpeed = 1000000L;  // 1 Mbps
  const double packetLength = 1500L * 8;  // 1500 bytes
  const unsigned long ticks = 1000000L;

  numConnected = atol(argv[1]);
  arrivalRate = (double) atol(argv[2]);
  persistence = (double) atol(argv[3]);

  Simulation sim(numConnected, ticks, arrivalRate, packetLength, lanSpeed,
                 persistence);

  sim.start();

  cout << "Average packet delay: " << sim.getPacketDelay() << endl;
  cout << "Throughput: " << sim.getThroughput() << endl;
  sim.printMiscStats();

}
