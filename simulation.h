#ifndef _SIMULATION_H_
#define _SIMULATION_H_
#include "Debug.h"
#include <queue>
#include <map>

struct Packet {
 public:
  Packet(unsigned long id, unsigned long startTime);
 private:
  friend class Simulation;
  friend class Node;
  unsigned long m_id;
  unsigned long m_startTime;
  unsigned long m_endTime;
  unsigned int m_retry;
  unsigned int m_sourceId;
  unsigned int m_destId;
  unsigned long m_propEnd;  // propagate until tick
  unsigned long m_startTx;
  unsigned long m_start;
  unsigned long m_end;
};

class Node {
 public:
  enum State {
    Idle,
    Transmitting,
    Backoff
  };
  Node(unsigned int id, double A, double mediumSense);
 private:
  friend class Simulation;
  unsigned int m_id;
  unsigned long m_packetsRemaining;  // queue
  unsigned long m_nextPacketArrival;
  unsigned long m_wait;
  unsigned long m_backoff;  // backoff value
  unsigned long m_collisions;
  double m_lambda;  // arrival rate
  double m_mediumSense;
  bool m_mediumBusy;
  State m_state;
  Packet *m_currentPacket;
  std::queue<Packet*> m_queue;
  std::queue<Packet*> m_failQueue;
  unsigned long m_busyEnd;  // busy until
  std::map<unsigned long, Packet*> m_propagatingPackets;
  std::map<unsigned long, Packet*> m_passingPackets;

  std::map<unsigned long, unsigned long> *detectBusyMedium(unsigned long ticks);
  unsigned long nextInterval(void);
  void generateNextPacketArrival(unsigned long tick);

};

class Simulation {
 public:
  Simulation(unsigned int connectedNodes, unsigned long ticks, double lambda,
             double packetLength, double transmissionRate, double persistence);
  void start(void);
  double getPacketDelay() const;
  double getThroughput() const;
  void printMiscStats() const;

 private:
  static unsigned long UniqueCounter;
  unsigned int m_connectedNodes;  // number of connected nodes
  unsigned long m_ticks;
  double m_lambda;  // arrival rate
  double m_packetLength;
  double m_transmissionRate;  // lan speed
  double m_persistence;

  unsigned long m_packetsFail;
  unsigned long m_packetsSuccess;
  unsigned long m_packetsDelay;
  unsigned long m_collisions;
  unsigned long m_successTrans;
  unsigned long m_totalPackets;
  bool m_simFinished;

  double m_tp;
  double m_mediumSense;
  double m_dProp;
  double m_dTrans;
  double m_jammingSignal;

  Node **nodes;

  void iteratePassing(std::map<unsigned long, unsigned long> *overlap);
  void detectCollisions(Node &node, unsigned long ticks);
  void propagatePackets(Node &node, unsigned long ticks);
  Packet *processPacketArrival(Node &node, unsigned long tick);
};

#endif
