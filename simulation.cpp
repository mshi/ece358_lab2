#include <vector>
#include <limits>
#include <iostream>
#include <cstdlib>
#include "simulation.h"
#include <cmath>
#include "random.h"
#include "assert.h"
#define KMax 10 // i-saturation
#define TicksPerSecond 1e6

using namespace std;

Packet::Packet(unsigned long id, unsigned long startTime)
    : m_id(id),
      m_startTime(startTime) {
  m_endTime = 0;
  m_retry = 0;
  m_sourceId = numeric_limits<unsigned int>::max();
  m_destId = numeric_limits<unsigned int>::max();
  m_propEnd = numeric_limits<unsigned int>::max();
  m_startTx = 0;
}

Node::Node(unsigned int id, double A, double mediumSense)
    : m_id(id),
      m_packetsRemaining(0),
      m_nextPacketArrival(1),
      m_wait(0),
      m_backoff(0),
      m_collisions(0),
      m_lambda(A),
      m_state(Idle),
      m_mediumSense(mediumSense) {
  m_busyEnd = 0;
  m_mediumBusy = false;
  m_currentPacket = NULL;
}

Simulation::Simulation(unsigned int N, unsigned long ticks, double A, double L,
                       double W, double P)
    : m_connectedNodes(N),
      m_ticks(ticks),
      m_lambda(A),
      m_packetLength(L),
      m_transmissionRate(W),
      m_persistence(P),
      m_packetsFail(0),
      m_packetsSuccess(0),
      m_packetsDelay(0),
      m_simFinished(false),
      m_collisions(0),
      m_successTrans(0),
      m_totalPackets(0) {

  const double S = 2e8;  // speed of light in copper
  const double D = 10;  // 10m distance between neighbours
  m_tp = (512 * TicksPerSecond) / W;
  m_mediumSense = (96 * TicksPerSecond) / W;
  m_dProp = TicksPerSecond * (D / S);
  m_dTrans = TicksPerSecond * (L / W);
  m_jammingSignal = 48 * TicksPerSecond;
  nodes = new Node*[N];
  cout << "dProp: " << m_dProp << endl;
  cout << "dTrans: " << m_dTrans << endl;
  cout << "Tp: " << m_tp << endl;
}

void Node::generateNextPacketArrival(unsigned long tick) {
  if (tick > m_nextPacketArrival) {
    // Generate time for next packet arrival
    m_nextPacketArrival = tick + nextInterval();
  }
}

unsigned long Node::nextInterval(void) {
  double u = genrand();
  return ceil(-1 / m_lambda * TicksPerSecond * log(u));
}

Packet *Simulation::processPacketArrival(Node &node, unsigned long tick) {
  if (node.m_nextPacketArrival == tick) {  // generate new packet
    Debug::getInstance() << "Generate new packet" << endl;
    ++m_totalPackets;
    Packet *p = new Packet(tick, tick);
    // randomly assign packet to stations
    unsigned int src = rand() % m_connectedNodes;  // (0,N)
    unsigned int dest = rand() % m_connectedNodes;
    while (src == dest) {
      dest = rand() % m_connectedNodes;
    }
    // put it in queue of the random station
    p->m_sourceId = nodes[src]->m_id;
    p->m_destId = dest;
    nodes[src]->m_queue.push(p);
    return p;
  }
  return NULL;
}

void Simulation::start(void) {

  Debug::getInstance() << "Start simulation" << endl;

  queue<Packet*> listOfPackets;

  for (unsigned int i = 0; i < m_connectedNodes; ++i) {
    Debug::getInstance() << "Creating node: " << i << endl;
    nodes[i] = new Node(i, m_lambda, m_mediumSense);
  }

  Packet *p;

  for (unsigned long currentTick = 1; currentTick < m_ticks; ++currentTick) {

    cout << "current tick: " << currentTick << endl;
    // go through all the nodes
    for (int i = 0; i < m_connectedNodes; ++i) {
      Node &currentNode = *(nodes[i]);

      currentNode.generateNextPacketArrival(currentTick);
      p = processPacketArrival(currentNode, currentTick);
      if (p != NULL) {
        listOfPackets.push(p);
      }

      switch (currentNode.m_state) {
        case Node::Transmitting:
          Debug::getInstance() << "Node " << currentNode.m_id
                               << " is transmitting" << endl;
          if (currentNode.m_busyEnd == currentTick) {  // if it finished transmitting
            currentNode.m_mediumBusy = false;
            ++m_successTrans;
            currentNode.m_state = Node::Idle;
            currentNode.m_mediumSense = m_mediumSense;
          } else {  // process and wait
            if (!currentNode.m_mediumBusy) {
              if (currentNode.m_failQueue.empty()) {
                currentNode.m_currentPacket = currentNode.m_queue.front();
                currentNode.m_queue.pop();
              } else {
                currentNode.m_currentPacket = currentNode.m_failQueue.front();
                currentNode.m_failQueue.pop();
              }
              currentNode.m_propagatingPackets[currentNode.m_currentPacket->m_id] =
                  currentNode.m_currentPacket;
              currentNode.m_currentPacket->m_startTx = currentTick;

              // compute passing packets
              int direction =
                  (currentNode.m_id > currentNode.m_currentPacket->m_destId) ?
                      -1 : 1;
              for (int nId = currentNode.m_id + direction;
                  nId != (currentNode.m_currentPacket->m_destId + direction);
                  nId = nId + direction) {
                unsigned long start = currentTick
                    + (abs(nId - currentNode.m_id)) * m_dProp;
                Packet *infoPacket = new Packet(0, start);
                infoPacket->m_endTime = start + m_dTrans;
                infoPacket->m_sourceId = currentNode.m_id;
                infoPacket->m_destId = currentNode.m_currentPacket->m_destId;
                nodes[nId]->m_passingPackets[currentNode.m_currentPacket->m_id] =
                    infoPacket;
                listOfPackets.push(infoPacket);
              }

              currentNode.m_mediumBusy = true;
              currentNode.m_busyEnd = currentTick + m_dTrans;
              currentNode.m_currentPacket->m_propEnd = currentTick
                  + (abs(
                      currentNode.m_currentPacket->m_destId - currentNode.m_id))
                      * m_dProp + m_dTrans;
              Debug::getInstance() << "Node: " << currentNode.m_id
                                   << " busy until " << currentNode.m_busyEnd
                                   << endl;
              Debug::getInstance() << "Node: " << currentNode.m_id
                                   << " propogating packet: "
                                   << currentNode.m_currentPacket->m_id
                                   << " until "
                                   << currentNode.m_currentPacket->m_propEnd
                                   << endl;
            }
          }
          detectCollisions(currentNode, currentTick);
          break;
        case Node::Idle:
//          Debug::getInstance() << "Node " << currentNode.m_id << " is idle"
//                               << endl;
          iteratePassing(currentNode.detectBusyMedium(currentTick));
          currentNode.m_wait = 0;
          if (currentNode.m_mediumBusy) {
            currentNode.m_wait = 10;
            currentNode.m_mediumSense = m_mediumSense;
          }
//          Debug::getInstance() << "Node " << currentNode.m_id << " is busy? "
//                               << currentNode.m_mediumBusy << " sense: "
//                               << currentNode.m_mediumSense << " wait: "
//                               << currentNode.m_wait << endl;
          if (currentNode.m_wait > 0) {
            --currentNode.m_wait;
          } else if (!currentNode.m_queue.empty()) {
            if (!currentNode.m_mediumBusy && currentNode.m_mediumSense <= 0) {
              Debug::getInstance() << "Changing to transmit for Node: "
                                   << currentNode.m_id << endl;
              currentNode.m_state = Node::Transmitting;
            }
            --currentNode.m_mediumSense;
          }
          break;
        case Node::Backoff:
          Debug::getInstance() << "Node " << currentNode.m_id
                               << " is in backoff wait" << endl;
          iteratePassing(currentNode.detectBusyMedium(currentTick));
          if (!currentNode.m_mediumBusy) {
            Debug::getInstance() << "Node: " << currentNode.m_id << " Backoff: "
                                 << currentNode.m_backoff << endl;
            --(currentNode.m_backoff);
            if (currentNode.m_backoff <= 0) {
              currentNode.m_mediumSense = m_mediumSense;
              currentNode.m_state = Node::Idle;
            }
          }
          break;
        default:
          Debug::getInstance() << "Bad news bears for Node: "
                               << currentNode.m_id << endl;
          break;
      }
      propagatePackets(currentNode, currentTick);
    }
  }

  m_simFinished = true;
  while (!listOfPackets.empty()) {
    Packet *toDelete = listOfPackets.front();
    listOfPackets.pop();
    delete toDelete;
  }
  for (unsigned int i = 0; i < m_connectedNodes; ++i) {
    delete nodes[i];
  }
  delete[] nodes;
}

void Simulation::iteratePassing(map<unsigned long, unsigned long> *overlap) {
  if (overlap == NULL) {
    return;
  }
  // iterate through colliding packets
  if (overlap->size() > 1) {
    for (map<unsigned long, unsigned long>::iterator src =
        overlap->begin(); src != overlap->end(); ++src) {
      unsigned long nodeSrc = src->first;
      ++(nodes[nodeSrc]->m_propagatingPackets[overlap->at(nodeSrc)]->m_retry);
      if (nodes[nodeSrc]->m_propagatingPackets[overlap->at(nodeSrc)]->m_retry
          > KMax) {
        ++m_packetsFail;
        nodes[nodeSrc]->m_state = Node::Idle;
        nodes[nodeSrc]->m_mediumBusy = false;
      } else {
        nodes[nodeSrc]->m_failQueue.push(
            nodes[nodeSrc]->m_propagatingPackets[overlap->at(nodeSrc)]);
        double r = genrand()
            * (pow(
                2,
                nodes[nodeSrc]->m_propagatingPackets[overlap->at(nodeSrc)]
                    ->m_retry) - 1);
        nodes[nodeSrc]->m_state = Node::Backoff;
        nodes[nodeSrc]->m_mediumBusy = false;
        nodes[nodeSrc]->m_backoff = r * m_tp;
      }
      // remove passing packets
      for (int i = 0; i < m_connectedNodes; ++i) {
        map<unsigned long, Packet*>::iterator it = nodes[i]->m_passingPackets
            .find(overlap->at(nodeSrc));
        if (it != nodes[i]->m_passingPackets.end())
          nodes[i]->m_passingPackets.erase(it);
      }
      map<unsigned long, Packet*>::iterator it = nodes[nodeSrc]
          ->m_propagatingPackets.find(overlap->at(nodeSrc));
      if (it != nodes[nodeSrc]->m_passingPackets.end()) {
        nodes[nodeSrc]->m_passingPackets.erase(it);
      }
    }
  }
  delete overlap;
}

void Simulation::detectCollisions(Node &node, unsigned long ticks) {
  if (node.m_passingPackets.empty()) {
    return;
  }

  for (map<unsigned long, Packet*>::iterator i = node.m_passingPackets
      .begin(); i != node.m_passingPackets.end(); ++i) {

    if(node.m_passingPackets.empty()){
      break;
    }
    if(i->second == NULL){
      break;
    }

    if ((ticks >= i->second->m_startTime) && (ticks < i->second->m_endTime)) {
      // packet collision retransmit packet
      ++m_collisions;
      ++(node.m_currentPacket->m_retry);
      if (node.m_currentPacket->m_retry > KMax) {
        ++m_packetsFail;
        node.m_state = Node::Idle;
        node.m_mediumBusy = false;
      } else {
        node.m_failQueue.push(node.m_currentPacket);
        double r = genrand() * (pow(2, node.m_currentPacket->m_retry) - 1);
        node.m_state = Node::Backoff;
        node.m_mediumBusy = false;
        node.m_backoff = r * m_tp;
      }
      map<unsigned long, Packet*>::iterator it = node.m_propagatingPackets.find(
          node.m_currentPacket->m_id);
      if (it != node.m_propagatingPackets.end()) {
        node.m_propagatingPackets.erase(it);
      }
      unsigned long collided = numeric_limits<unsigned long>::max();
      // retransmit packet that caused collision
      unsigned long nodeSrc = i->second->m_sourceId;
      for (map<unsigned long, Packet*>::iterator j = nodes[nodeSrc]
          ->m_propagatingPackets.begin();
          j != nodes[nodeSrc]->m_propagatingPackets.end(); ++j) {
        if (j->first == i->first) {
          ++(j->second->m_retry);
          if (j->second->m_retry > KMax) {
            ++m_packetsFail;
            nodes[nodeSrc]->m_state = Node::Idle;
            nodes[nodeSrc]->m_mediumBusy = false;
          } else {
            nodes[nodeSrc]->m_failQueue.push(j->second);
            double r = genrand()
                * (pow(2, nodes[nodeSrc]->m_currentPacket->m_retry) - 1);
            nodes[nodeSrc]->m_state = Node::Backoff;
            nodes[nodeSrc]->m_mediumBusy = false;
            nodes[nodeSrc]->m_backoff = r * m_tp;
          }
          collided = j->first;
        }
      }
      // delete propogating packets
      if (collided != numeric_limits<unsigned long>::max()) {
        it = nodes[nodeSrc]->m_propagatingPackets.find(collided);
        if (it != nodes[nodeSrc]->m_propagatingPackets.end()) {
          nodes[nodeSrc]->m_propagatingPackets.erase(it);
        }
      }
      for (int n = 0; n < m_connectedNodes; ++n) {
        // packet that caused collision
        it = nodes[n]->m_passingPackets.find(collided);
        if (it != nodes[n]->m_passingPackets.end()) {
          nodes[n]->m_passingPackets.erase(it);
        }

        // packet that was transmitting
        it = nodes[n]->m_passingPackets.find(node.m_currentPacket->m_id);
        if (it != nodes[n]->m_passingPackets.end()) {
          nodes[n]->m_passingPackets.erase(it);
        }
      }
    }
  }
}

map<unsigned long, unsigned long> *Node::detectBusyMedium(unsigned long ticks) {

  if (m_passingPackets.empty()) {
    m_mediumBusy = false;
    return NULL;
  }

  map<unsigned long, unsigned long> *overlappingPackets = new map<unsigned long,
      unsigned long>();
  for (map<unsigned long, Packet*>::const_iterator i = m_passingPackets.begin();
      i != m_passingPackets.end(); ++i) {
    if ((ticks >= i->second->m_startTime) && (ticks < i->second->m_endTime)) {  // overlapping
      overlappingPackets->insert(
          make_pair<unsigned long, unsigned long>(i->second->m_sourceId,
                                                  i->first));
      m_mediumBusy = true;
    } else {  // not overlapping
      m_mediumBusy = false;
    }
  }

  return overlappingPackets;
}

void Simulation::propagatePackets(Node &node, unsigned long ticks) {
  vector<unsigned long> received;

  for (map<unsigned long, Packet*>::iterator i =
      node.m_propagatingPackets.begin(); i != node.m_propagatingPackets.end();
      ++i) {
    if (ticks >= i->second->m_propEnd) {
      received.push_back(i->first);
      node.m_mediumBusy = false;
      i->second->m_endTime = ticks;
      ++m_packetsSuccess;
      m_packetsDelay += ticks - i->second->m_startTime;
    }
  }
  if (!received.empty()) {
    for (vector<unsigned long>::const_iterator i = received.begin();
        i != received.end(); ++i) {
      map<unsigned long, Packet*>::iterator it = node.m_propagatingPackets.find(
          *i);
      if (it != node.m_propagatingPackets.end()) {
        node.m_propagatingPackets.erase(it);
      }
      for (int j = 0; j < m_connectedNodes; ++j) {
        it = nodes[j]->m_passingPackets.find(*i);
        if (it != nodes[j]->m_passingPackets.end()) {
          nodes[j]->m_passingPackets.erase(it);
        }
      }
    }
  }
}

double Simulation::getPacketDelay() const {
  return m_packetsDelay / m_packetsSuccess;
}

double Simulation::getThroughput() const {
  return (m_packetsSuccess * m_packetLength * TicksPerSecond)
      / (m_ticks * m_transmissionRate);
}

void Simulation::printMiscStats() const {
  cout << "Successful receives: " << m_packetsSuccess << endl;
  cout << "Successful transmit: " << m_successTrans << endl;
  cout << "Failed packets: " << m_packetsFail << endl;
  cout << "Collided Packets: " << m_collisions << endl;
  cout << "Total Packets: " << m_totalPackets << endl;

}
