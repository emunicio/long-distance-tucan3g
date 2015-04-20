/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *      Author: Esteban Municio <esteban.municio@urjc.es>
 *
 ******************************************************************/

/*
This script simulate a WiMAX Long Distance link in a simple manner. It was developed for obtain throughput and latency values when the link is working in the staturation point for different distances with a specific link configuration (Modulation, Packet Size, Frame Duration, Cyclic Prefix, Bandwidth)

      +-----+			+------------+
      | SS0 |   (<---------)    |Base Station|
      +-----+          	        +------------+
      10.1.1.1                    10.1.1.2

The accepted paramterers are:

	--mod:       Modulation and codification scheme used. Different modulation schemes can be used, 
		     from BPSK1/2 to 64QAM3/4. Only SISO techniques are considered.
	--frame:     Frame Duration used in seconds. The following frame durations are allowed:
		     0.0025, 0.004, 0.005, 0.008, 0.010, 0.012 and 0.020 seconds.
	--cp:	     Cyclic Prefix used. The following Cyclic Prefix are allowed:
		     0.24, 0.125, 0.0625 and 0.03125.
	--pSize:     Fixed packet size used for the simulation in bytes. This value correspond to APP 
		     level size. It is necessary to consider the protocol overhead since results are given at 
		     MAC level.
	--bw:	     Bandwidth used in Hz. Only the standard [4] bandwiths for are allowed: 7 and 10 MH>
		     for licensed bands and non-licensed bands respectively.
	--pcap:      Generate Pcap trace file. If true, generate a Pcap file.

The following default values are used in this scripts:
	Standard:		IEEE 802.16e-2007
	Frequency:		5 GHz
	Distance:		8 Km
	UL/DL:			0%
	Scheduler:		Simple
	RTG/TTG:		Distance dependant
	Ranging:		Distance dependant
	WimaxMacQueue		1024
	QoS:			Service Flow UGS
	Transport Protocol:	UDP
	Propagation Model:	Random Propagation
	Bidirectional flows:	No

Example of use:

./waf --run "scratch/basic-WiMAX --mod=0 --frame=0.02 --cp=0.25 --pSize=1300 --bw=10000000 --pcap=true"

In this experiment, results of throughput and latency are given for the saturation state. The distance used is fixed to 8 Km. The modulation used is BPSK 1/2 (0), the frame duration is 20ms, the Cyclic Prefix is 1/4 and the packet size at APP level is 1300 bytes. The bandwidth considered is 10 MHz. 2 pcap files will be generated for both BS and SS.
*/



#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/applications-module.h"
#include "ns3/mobility-module.h"
#include "ns3/config-store-module.h"
#include "ns3/wimax-module.h"
#include "ns3/internet-module.h"
#include "ns3/ipcs-classifier-record.h"
#include "ns3/service-flow.h"
#include <iostream>
#include "ns3/mobility-module.h"
#include "ns3/global-route-manager.h"
#include "ns3/service-flow.h"
#include "ns3/flow-monitor-module.h"

using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("TUCAN3G");

class Experiment
{
	public:
		Experiment ();
		void RunSim (uint32_t dataRate, uint32_t distance, uint32_t packetSize, uint32_t modulation,
				  double frame_duration, double cp, uint32_t bandwidth, uint32_t runApplicationSeconds,
					int32_t rssi, bool pcapTrace);
		uint64_t GetThroughput ();
		uint64_t GetDelay ();

		void PhyRxDl (std::string context , ns3::Ptr<ns3::Packet const> p,ns3::Mac48Address const &mac);
		void PhyRxUl (std::string context , ns3::Ptr<ns3::Packet const> p,ns3::Mac48Address const &mac);
		void PhyTxDl (std::string context , ns3::Ptr<ns3::Packet const> p,ns3::Mac48Address const &mac);
		void PhyTxUl (std::string context , ns3::Ptr<ns3::Packet const> p,ns3::Mac48Address const &mac);


	private:

		Ptr<Socket> SetupPacketReceive1 (Ipv4Address addr, Ptr<Node> node);
		Ptr<Socket> SetupPacketReceive2 (Ipv4Address addr, Ptr<Node> node);

		void ReceivePacket1 (Ptr<Socket> socket);
		void ReceivePacket2 (Ptr<Socket> socket);

		uint32_t port1;
		uint32_t port2;
		uint32_t m_dataRate;
		uint32_t m_distance;
		uint32_t m_packetSize;


		float m_delay;
		uint32_t m_runApplicationSeconds;
		int32_t m_rssi;
		uint32_t m_modulation;
		double m_frame_duration;
		double m_cp;

		uint32_t m_RxPacketsDl;
		uint32_t m_TxPacketsDl;
		uint32_t m_RxPacketsUl;
		uint32_t m_TxPacketsUl;
		uint32_t m_RxBytesDl;
		uint32_t m_TxBytesDl;
		uint32_t m_RxBytesUl;
		uint32_t m_TxBytesUl;
		double m_delayDl;
		double m_delayUl;
		uint32_t m_throughput;
		uint32_t m_bandwidth;

};


Experiment::Experiment () :
		port1(4002),
		port2(4000),
		m_RxPacketsDl(0),
		m_TxPacketsDl(0),
		m_RxPacketsUl(0),
		m_TxPacketsUl(0),
		m_RxBytesDl(0),
		m_TxBytesDl(0),
		m_RxBytesUl(0),
		m_TxBytesUl(0),
		m_delayDl(0),
		m_delayUl(0),
		m_throughput(0),
		m_bandwidth(0)
{	}

void
Experiment::ReceivePacket1 (Ptr<Socket> socket)
{

  Ptr<Packet> packet;
    while (packet = socket->Recv ())
    {
      //bytesTotal1 += (packet->GetSize()+8+20+8);
      //packetsReceived1 += 1;
    }
 }

void
Experiment::ReceivePacket2 (Ptr<Socket> socket)
{

  Ptr<Packet> packet;
    while (packet = socket->Recv ())
    {
      //bytesTotal2 += (packet->GetSize()+8+20+8);
      //packetsReceived2 += 1;
    }
 }

Ptr<Socket>
Experiment::SetupPacketReceive1 (Ipv4Address addr, Ptr<Node> node)
{

  TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
  Ptr<Socket> sink = Socket::CreateSocket (node, tid);
  InetSocketAddress local = InetSocketAddress (addr, port1);
  sink->Bind (local);
  sink->SetRecvCallback (MakeCallback (&Experiment::ReceivePacket1, this));
  return sink;
}

Ptr<Socket>
Experiment::SetupPacketReceive2 (Ipv4Address addr, Ptr<Node> node)
{

  TypeId tid = TypeId::LookupByName ("ns3::UdpSocketFactory");
  Ptr<Socket> sink = Socket::CreateSocket (node, tid);
  InetSocketAddress local = InetSocketAddress (addr, port2);
  sink->Bind (local);
  sink->SetRecvCallback (MakeCallback (&Experiment::ReceivePacket2, this));
  return sink;
}

void
  Experiment::PhyTxDl (std::string context , ns3::Ptr<ns3::Packet const> p,ns3::Mac48Address const &mac)
  {
	UdpHeader udph;
	p->PeekHeader(udph);
	//std::cout << "Time "<< Simulator::Now().GetSeconds () << " Size" << p->GetSize() << " DL PHY TX " << p->GetReferenceCount() << " " << p->GetUid() << std::endl;
	m_TxPacketsDl=m_TxPacketsDl+1;
	m_TxBytesDl=m_TxBytesDl+p->GetSize();

  }

void
  Experiment::PhyTxUl (std::string context , ns3::Ptr<ns3::Packet const> p,ns3::Mac48Address const &mac)
  {
	MacHeaderType hdr1;
	p->PeekHeader(hdr1);
	//std::cout << p->GetSize() << " UL PHY TX 1 from " << mac << std::endl;
	m_TxPacketsUl=m_TxPacketsUl+1;
	m_TxBytesUl=m_TxBytesUl+p->GetSize();

  }

  void
    Experiment::PhyRxDl (std::string context , ns3::Ptr<ns3::Packet const> p,ns3::Mac48Address const &mac)
    {
	   MacHeaderType hdr1;
	   p->PeekHeader(hdr1);
	   m_RxPacketsDl=m_RxPacketsDl+1;
	   m_RxBytesDl=m_RxBytesDl+p->GetSize();
	   //std::cout <<"["<<Simulator::Now ().GetSeconds()<<"]"<< " RX time" <<  std::endl;
	   TimeStampOnOff timestamp;
	   p->FindFirstMatchingByteTag(timestamp);
	   //std::cout <<"["<< timestamp.GetTimestamp().GetSeconds() <<"]"<< " Tx time " << timestamp.GetTypeId() << std::endl;
	   m_delayDl=m_delayDl+Simulator::Now ().GetSeconds()-timestamp.GetTimestamp().GetSeconds ();
	   //std::cout <<m_delayDl<<  std::endl;
    }
  void
      Experiment::PhyRxUl (std::string context , ns3::Ptr<ns3::Packet const> p,ns3::Mac48Address const &mac)
      {
		MacHeaderType hdr1;
		p->PeekHeader(hdr1);
		m_RxPacketsUl=m_RxPacketsUl+1;
		m_RxBytesUl=m_RxBytesUl+p->GetSize();

		TimeStampOnOff timestamp;
		p->FindFirstMatchingByteTag(timestamp);
		//std::cout << "Rx packet with delay "<< Simulator::Now ().GetSeconds()-timestamp.GetTimestamp().GetSeconds ()<< " Size:  "<< p->GetSize() << " UL PHY RX " << p->GetReferenceCount() << " "<< p->GetUid() << std::endl;
		//m_delayUl=m_delayUl+Simulator::Now ().GetSeconds()-timestamp.GetTimestamp().GetSeconds ();
      }



void
PopulateArpCache ()
{
  Ptr<ArpCache> arp = CreateObject<ArpCache> ();
  arp->SetAliveTimeout (Seconds(3600 * 24 * 365));
  for (NodeList::Iterator i = NodeList::Begin(); i != NodeList::End(); ++i)
  {
    Ptr<Ipv4L3Protocol> ip = (*i)->GetObject<Ipv4L3Protocol> ();
    NS_ASSERT(ip !=0);
    ObjectVectorValue interfaces;
    ip->GetAttribute("InterfaceList", interfaces);
    for(ObjectVectorValue::Iterator j = interfaces.Begin(); j !=
interfaces.End (); j ++)
    {
    	Ptr<Ipv4Interface> ipIface = (j->second)->GetObject<Ipv4Interface> ();
      NS_ASSERT(ipIface != 0);
      Ptr<NetDevice> device = ipIface->GetDevice();
      NS_ASSERT(device != 0);
      Mac48Address addr = Mac48Address::ConvertFrom(device->GetAddress ());
      for(uint32_t k = 0; k < ipIface->GetNAddresses (); k ++)
      {
        Ipv4Address ipAddr = ipIface->GetAddress (k).GetLocal();
        if(ipAddr == Ipv4Address::GetLoopback())
          continue;
        ArpCache::Entry * entry = arp->Add(ipAddr);
        entry->MarkWaitReply(0);
        entry->MarkAlive(addr);
        //NS_LOG_UNCOND ("Arp Cache: Adding the pair (" << addr << "," << ipAddr << ")");
      }
    }
  }
  for (NodeList::Iterator i = NodeList::Begin(); i != NodeList::End(); ++i)
  {
    Ptr<Ipv4L3Protocol> ip = (*i)->GetObject<Ipv4L3Protocol> ();
    NS_ASSERT(ip !=0);
    ObjectVectorValue interfaces;
    ip->GetAttribute("InterfaceList", interfaces);
    for(ObjectVectorValue::Iterator j = interfaces.Begin(); j !=interfaces.End (); j ++)
    {
    	Ptr<Ipv4Interface> ipIface = (j->second)->GetObject<Ipv4Interface> ();
      ipIface->SetAttribute("ArpCache", PointerValue(arp));
    }
  }
}

void
PrintHelp ()
{
	 std::cout <<"\nProgram Arguments: \n";
	 std::cout <<"   --mod:       Modulation and codification scheme used [0]\n";
	 std::cout <<"   --frame:     Frame duration used in seconds [0.0025]\n";
	 std::cout <<"   --cp:        Cyclic Prefix used [0.25]\n";
	 std::cout <<"   --pSize:     Fixed packet size used for the simulation in bytes [1372]\n";
	 std::cout <<"   --bw:   	  Bandwidth used in Hz [10000000]\n";
	 std::cout <<"   --pcap:      Generate Pcap trace file (boolean) [false]\n";
}

  int main (int argc, char *argv[])
  {

  	Experiment experiment;

  	// Command line arguments

  	int32_t modulation = 0;	//mod 0-6

  	double frameDuration=0.0025;	//aggregation threshold
  	double cp=0.25;	//adjust ackTimeout for optimal behavior

  	// default values for type conversion
    std::string s_frameDuration="0.0025";	//in string
  	std::string s_cp="0.25";				//in string
  	int32_t i_frameDuration=9;				//in index integer
  	int32_t i_cp=9;							//in index integer

  	int32_t distance=8000; //max distance in meters
  	int32_t packetSize = 1372;	//Fixed packet size used for the simulation (app level)
  	int32_t bandwidth=10000000;
    bool pcapTrace = false;
    int32_t rssi = 90; //received signal

  	CommandLine cmd;
  	cmd.AddValue ("mod", "Modulation and codification scheme used", modulation);
  	cmd.AddValue ("frame", "Frame Duration", s_frameDuration);
  	cmd.AddValue ("cp", "Cyclic prefix", s_cp);
  	cmd.AddValue ("pSize", "Fixed packet size used for the simulation in bytes", packetSize);
  	cmd.AddValue ("bw", "Bandwidth used ", bandwidth);
    cmd.AddValue ("pcap", "Create Pcap trace file (boolean)", pcapTrace);
  	cmd.Parse (argc, argv);

    std::cout << "\n******************************************************************";
  	std::cout << "\nBeging simulation with the following configuration:";

  	std::cout << "\nWiMAX Link: 8 Km";

  	switch (modulation)
  		     {
  		     case 0:
  		    	std::cout << "\nBPSK 1/2 -> " << modulation << " (Using Default)";
  		       break;
  		     case 1:
  		    	std::cout << "\nQPSK 1/2 -> " << modulation;
  		       break;
  		     case 2:
  		    	std::cout << "\nQPSK 3/4 -> " << modulation;
  		       break;
  		     case 3:
  		        std::cout << "\n16QAM 1/2 -> " << modulation;
  		       break;
  		     case 4:
  		   		std::cout << "\n16QAM 3/4 -> " << modulation;
  		   	   break;
  		     case 5:
  		  		std::cout << "\n64QAM 2/3 -> " << modulation;
  		  	   break;
  		     case 6:
  		     	std::cout << "\n64QAM 3/4 -> " << modulation;
  		       break;
  		     default:
  		    	std::cout << "\nBPSK 1/2 -> " << modulation;
				std::cout << "\nQPSK 1/2 -> " << modulation;
				std::cout << "\nQPSK 3/4 -> " << modulation;
				std::cout << "\n16QAM 1/2 -> " << modulation;
				std::cout << "\n16QAM 3/4 -> " << modulation;
				std::cout << "\n64QAM 2/3 -> " << modulation;
				std::cout << "\n64QAM 3/4 -> " << modulation;
				std::cout << "\nUsing:";
  		    	std::cout << "\nBPSK 1/2 -> " << modulation << " (Using Default)";
  		    	modulation=0;
  		     }

  	if(s_frameDuration == "0.0025"){
  		i_frameDuration=0;
  	}if(s_frameDuration == "0.004"){
  		i_frameDuration=1;
	}if(s_frameDuration == "0.005"){
		i_frameDuration=2;
	}if(s_frameDuration == "0.008"){
		i_frameDuration=3;
	}if((s_frameDuration == "0.010") || (s_frameDuration == "0.01")){
		i_frameDuration=4;
	}if((s_frameDuration == "0.012")|| (s_frameDuration == "0.0125")){
		i_frameDuration=5;
	}if((s_frameDuration == "0.02") || (s_frameDuration == "0.020")){
		i_frameDuration=6;
	}

	switch(i_frameDuration)
    {
		 case 0:
			 std::cout << "\nFrame Duration -> 0.0025 s (Using Default)";
			 frameDuration=0.0025;
		   break;
		 case 1:
			 std::cout << "\nFrame Duration -> 0.004 s";
			 frameDuration=0.004;
		   break;
		 case 2:
			 std::cout << "\nFrame Duration -> 0.005 s";
			 frameDuration=0.005;
		   break;
		 case 3:
			 std::cout << "\nFrame Duration -> 0.008 s";
			 frameDuration=0.008;
		   break;
		 case 4:
			 std::cout << "\nFrame Duration -> 0.010 s";
			 frameDuration=0.01;
		   break;
		 case 5:
			 std::cout << "\nFrame Duration -> 0.012 s";
			 frameDuration=0.0125;
		   break;
		 case 6:
			 std::cout << "\nFrame Duration -> 0.020 s";
			 frameDuration=0.02;
		   break;
		 default:
			std::cout << "\nOnly the following Frame Duration are allowed: -> ";
			std::cout << "\n0.0025, 0.004, 0.005, 0.008, 0.010, 0.012 and 0.020";
			std::cout << "\nUsing:";
			std::cout << "\nFrame Duration -> 0.0025 (Default)";
			frameDuration=0.0025;
     }


  			 if(s_cp == "0.25"){
  			  	i_cp=0;
  			 }if(s_cp == "0.125"){
  				i_cp=1;
			 }if(s_cp == "0.0625"){
				i_cp=2;
			 }if(s_cp == "0.03125"){
				i_cp=3;
			 }

  			switch (i_cp)
			 {
			 case 0: // 1/4
				 std::cout << "\nCyclic Prefix -> 1/4  (Using Default)";
				 cp=0.25;
			   break;
			 case 1: // 1/8
				 std::cout << "\nCyclic Prefix -> 1/8";
				 cp=0.125;
			   break;
			 case 2: // 1/16
				 std::cout << "\nCyclic Prefix -> 1/16";
				 cp=0.0625;
			   break;
			 case 3: // 1/32
				 std::cout << "\nCyclic Prefix -> 1/32";
				 cp=0.03125;
			   break;
			 default:
				std::cout << "\nOnly the following Cyclic Prefix values are allowed: -> ";
				std::cout << "\n0.25, 0.125, 0.0625 and 0.03125";
				std::cout << "\n Using:";
				std::cout << "\nCyclic Prefix -> 1/4 (Default)";
				cp=0.25;
			}

        if(packetSize==1372){
    	        std::cout << "\nPacketSize -> " << packetSize << " (Using Default)";
        }else{
    	        std::cout << "\nPacketSize -> " << packetSize;
        }

        switch (bandwidth)
		 {
		 case 7000000: // only for licensed bands
			 std::cout << "\nBandwidth -> 7 MHz";
		   break;
		 case 10000000: // only for non-licensed bands
			 std::cout << "\nBandwidth -> 10 MHz  (Using Default)";
		   break;
		 default:
			std::cout << "\nOnly the following Bandwidth values are allowed: -> ";
			std::cout << "\n 7000000 and 10000000";
			std::cout << "\n Using:";
			std::cout << "\nBandwidth -> 10 MHz (Default)";
			bandwidth=10000000;
		  }


        std::cout << "\n******************************************************************\n";

  	if(packetSize <= 0){
  	  	NS_LOG_UNCOND ("The packet size must be positive");
  	    PrintHelp();
  	  	return 0;
  	}

  	//simulation parameters

  	uint32_t runApplicationSeconds = 15;


  	// injected datarates
  	uint32_t dataRateL[7]={5200000,10000000,15000000,20000000,300000000,40000000,46000000};


    std::stringstream aa;
	aa << modulation << "";
	std::string mod = aa.str();

	std::string file_delay = "wimax_mod" + mod + "_frame"+ s_frameDuration +"_delay.ods";
	std::string file_throughput = "wimax_mod" + mod + "_frame"+ s_frameDuration +"_throughput.ods";



    //preparing output files
  	std::ofstream out_d1(file_delay.c_str(), std::ios::app);
  	out_d1 << "Distance (m) , Delay (ms)"<<	std::ends;
  	out_d1 << std::endl;
  	out_d1.close();

  	std::ofstream out_th1(file_throughput.c_str(), std::ios::app);
  	out_th1 << "Distance (m) , Sat. Throughput (bps)"<<	std::ends;
  	out_th1 << std::endl;
  	out_th1.close();

	experiment = Experiment();
	experiment.RunSim(dataRateL[modulation], distance, packetSize, modulation, frameDuration,cp, bandwidth, runApplicationSeconds,rssi,  pcapTrace);


   //fill output files with output values
   std::stringstream ss;
   ss <<distance;
   std::string sdistance = ss.str();
   std::ofstream out_d2(file_delay.c_str(), std::ios::app);
   out_d2 << sdistance <<" ,"<< experiment.GetThroughput() << " ,"<<	std::ends;
   out_d2 << std::endl;
   out_d2.close();

   std::ofstream out_th2(file_throughput.c_str(), std::ios::app);
   out_th2 << sdistance <<" ,"<< experiment.GetThroughput() << " ,"<<	std::ends;
   out_th2 << std::endl;
   out_th2.close();

	std::cout << "\n ---------------------------------------------";
	std::cout << "\nDelay values are in: \n" << file_delay;
	std::cout << "\nThroughput values are in: \n" << file_throughput;
	if(pcapTrace){
		std::cout << "\nBS Pcap traces are in: wimax-simple-bs0.pcap";
		std::cout << "\nSS Pcap traces are in: wimax-simple-ss0.pcap \n";
	}
    std::cout << "\nDone..." <<std::endl;
  	return 0;
  }

  uint64_t Experiment::GetThroughput (){
	 return m_throughput;
  }
  uint64_t Experiment::GetDelay (){
  	 return m_delay;
    }

 void
 Experiment::RunSim (uint32_t dataRate, uint32_t distance, uint32_t packetSize, uint32_t modulation,
		  double frame_duration, double cp, uint32_t bandwidth, uint32_t runApplicationSeconds,
			int32_t rssi, bool pcapTrace)
  {

	 	 	m_dataRate=dataRate;
	 	 	m_distance = distance;
	 	 	m_frame_duration=frame_duration;
	 	 	m_cp=cp;
	 	 	m_packetSize = packetSize;
	 	    m_modulation=modulation;
	 	 	m_runApplicationSeconds = runApplicationSeconds;
	 	 	m_bandwidth=bandwidth;


	 	 	int  schedType = 0;

	 	 	std::cout << "\nSimulating with: -------------------------------------------------------";
	 	 	std::cout << "\nModulation -> " << m_modulation;
	 	 	std::cout << "\nFrame Duration -> " << m_frame_duration;
	 	 	std::cout << "\nPacketSize -> " << m_packetSize;
	 	 	std::cout << "\nDistance -> " << m_distance;
	 	 	std::cout << "\nCP -> " << m_cp;
	 	 	std::cout << "\nBandwidth -> " << m_bandwidth;


	   WimaxHelper::SchedulerType scheduler = WimaxHelper::SCHED_TYPE_SIMPLE;
	   switch (schedType)
	     {
	     case 0:
	       scheduler = WimaxHelper::SCHED_TYPE_SIMPLE;
	       break;
	     case 1:
	       scheduler = WimaxHelper::SCHED_TYPE_MBQOS;
	       break;
	     case 2:
	       scheduler = WimaxHelper::SCHED_TYPE_RTPS;
	       break;
	     default:
	       scheduler = WimaxHelper::SCHED_TYPE_SIMPLE;
	     }

	   NodeContainer ssNode;
	   NodeContainer bsNode;

	   ssNode.Create (1);
	   bsNode.Create (1);

	   WimaxHelper wimax;


	   Config::SetDefault ("ns3::SimpleOfdmWimaxPhy::G", DoubleValue (m_cp));
	   Config::SetDefault ("ns3::SimpleOfdmWimaxPhy::TxPower", DoubleValue (rssi));
	   Config::SetDefault ("ns3::WimaxPhy::FrameDuration", TimeValue (Seconds (m_frame_duration)));
	   Config::SetDefault ("ns3::WimaxPhy::Bandwidth", UintegerValue (m_bandwidth));
	   Config::SetDefault ("ns3::WimaxMacQueue::MaxSize", UintegerValue(1024));


	   wimax.SetPropagationLossModel (SimpleOfdmWimaxChannel::RANDOM_PROPAGATION);


	   NetDeviceContainer ssDev, bsDev;

	   ssDev = wimax.Install (ssNode,
	                           WimaxHelper::DEVICE_TYPE_SUBSCRIBER_STATION,
	                           WimaxHelper::SIMPLE_PHY_TYPE_OFDM,
	                           scheduler);
	   bsDev = wimax.Install (bsNode, WimaxHelper::DEVICE_TYPE_BASE_STATION, WimaxHelper::SIMPLE_PHY_TYPE_OFDM, scheduler);


	   Ptr<SubscriberStationNetDevice> ss;

	   ss = ssDev.Get (0)->GetObject<SubscriberStationNetDevice> ();

	       switch (m_modulation)
	           {
	           case 0:
	         	  ss->SetModulationType (WimaxPhy::MODULATION_TYPE_BPSK_12);
	             break;
	           case 1:
	         	  ss->SetModulationType (WimaxPhy::MODULATION_TYPE_QPSK_12);
	             break;
	           case 2:
	         	  ss->SetModulationType (WimaxPhy::MODULATION_TYPE_QPSK_34);
	             break;
	           case 3:
	         	  ss->SetModulationType (WimaxPhy::MODULATION_TYPE_QAM16_12);
	             break;
	           case 4:
	         	  ss->SetModulationType (WimaxPhy::MODULATION_TYPE_QAM16_34);
	             break;
	           case 5:
	         	  ss->SetModulationType (WimaxPhy::MODULATION_TYPE_QAM64_23);
	             break;
	           case 6:
	         	  ss->SetModulationType (WimaxPhy::MODULATION_TYPE_QAM64_34);
	             break;
	           default:
	             NS_FATAL_ERROR ("Invalid modulation type");
	             break;
	           }

	   Ptr<BaseStationNetDevice> bs;

	   bs = bsDev.Get (0)->GetObject<BaseStationNetDevice> ();


		//mobility
	   Ptr<MobilityModel> bsPosition = CreateObject<ConstantPositionMobilityModel>();
	   bsPosition->SetPosition (Vector(0,0,0));
	   bsNode.Get(0)->AggregateObject (bsPosition);

	   Ptr<MobilityModel> SSPosition;
	   //Service stations
	   SSPosition = CreateObject<ConstantPositionMobilityModel> ();
	   SSPosition->SetPosition (Vector(m_distance,0,0));
	   ssNode.Get (0)-> AggregateObject (SSPosition);

	   InternetStackHelper stack;
	   stack.Install (bsNode);
	   stack.Install (ssNode);

	   Ipv4AddressHelper address;
	   address.SetBase ("10.1.1.0", "255.255.255.0");

	   Ipv4InterfaceContainer SSinterface = address.Assign (ssDev);
	   Ipv4InterfaceContainer BSinterface = address.Assign (bsDev);

	   if (pcapTrace)
	     {
		   wimax.EnablePcap ("wimax-simple-bs0", bsNode.Get (0)->GetId (), bs->GetIfIndex ());
		   wimax.EnablePcap ("wimax-simple-ss0", ssNode.Get (0)->GetId (), bs->GetIfIndex ());
		   //wimax.EnableLogComponents ();
	     }

	  OnOffHelper onoff2 ("ns3::UdpSocketFactory", InetSocketAddress (BSinterface.GetAddress (0), port1));
	  onoff2.SetAttribute ("OnTime", StringValue ("ns3::ConstantRandomVariable[Constant=1]"));
	  onoff2.SetAttribute ("OffTime", StringValue ("ns3::ConstantRandomVariable[Constant=0]"));
	  onoff2.SetAttribute ("DataRate", DataRateValue (DataRate (dataRate)));
	  onoff2.SetAttribute ("PacketSize", UintegerValue (m_packetSize));
	  //onoff2.SetAttribute("AccessClass",UintegerValue (0));
	  AddressValue remoteAddress2 (InetSocketAddress (SSinterface.GetAddress (0), port2));
	  onoff2.SetAttribute ("Remote", remoteAddress2);

	  ApplicationContainer temp2 = onoff2.Install (bsNode.Get(0));

	  temp2.Start(Seconds (2));//evito que haya colisiones al principio
      temp2.Stop (Seconds (2 + m_runApplicationSeconds));

	  Ptr<Socket> sink2 = SetupPacketReceive2 (SSinterface.GetAddress (0), ssNode.Get(0));


	 /*OnOffHelper onoff1 ("ns3::UdpSocketFactory", InetSocketAddress (SSinterface.GetAddress (0), port1));
	  onoff1.SetAttribute ("OnTime", StringValue ("ns3::ConstantRandomVariable[Constant=1]"));
	  onoff1.SetAttribute ("OffTime", StringValue ("ns3::ConstantRandomVariable[Constant=0]"));
	  onoff1.SetAttribute ("DataRate", DataRateValue (DataRate (m_datarate)));
	  onoff1.SetAttribute ("PacketSize", UintegerValue (m_packetSize));
	  //onoff2.SetAttribute("AccessClass",UintegerValue (0));
		AddressValue remoteAddress1 (InetSocketAddress (BSinterface.GetAddress (0), port2));
		onoff1.SetAttribute ("Remote", remoteAddress1);

		ApplicationContainer temp1 = onoff1.Install (ssNode.Get(0));

		temp1.Start(Seconds (warming));//evito que haya colisiones al principio
		temp1.Stop (Seconds (warming + m_runApplicationSeconds));


	  Ptr<Socket> sink1 = SetupPacketReceive1 (BSinterface.GetAddress (0), bsNode.Get(0));*/



	  Config::Connect ("/NodeList/0/DeviceList/0/$ns3::WimaxNetDevice/Rx",MakeCallback(&Experiment::PhyRxDl,this));
	  Config::Connect ("/NodeList/0/DeviceList/0/$ns3::WimaxNetDevice/Tx",MakeCallback(&Experiment::PhyTxUl,this));
	  Config::Connect ("/NodeList/1/DeviceList/0/$ns3::WimaxNetDevice/Rx",MakeCallback(&Experiment::PhyRxUl,this));
	  Config::Connect ("/NodeList/1/DeviceList/0/$ns3::WimaxNetDevice/Tx",MakeCallback(&Experiment::PhyTxDl,this));


	  Simulator::Stop (Seconds (2 + m_runApplicationSeconds+1));


	  //Flow configuration

	   IpcsClassifierRecord DlClassifierUgs (BSinterface.GetAddress (0),
	 		  	  	  	  	  	  	  	  Ipv4Mask ("255.255.255.255"),
	                                         SSinterface.GetAddress (0),
	                                         Ipv4Mask ("255.255.255.255"),
	                                         0,
	                                         65000,
	                                         port2,
	                                         port2,
	                                         17,
	                                         3);

	   ServiceFlow DlServiceFlowUgs = wimax.CreateServiceFlow (ServiceFlow::SF_DIRECTION_DOWN,
	                                                           ServiceFlow::SF_TYPE_UGS,
	                                                           DlClassifierUgs);

	   IpcsClassifierRecord UlClassifierUgs (SSinterface.GetAddress (0),
	                                         Ipv4Mask ("255.255.255.255"),
	                                         BSinterface.GetAddress (0),
	                                         Ipv4Mask ("255.255.255.255"),
	                                         0,
	                                         65000,
	                                         port1,
	                                         port1,
	                                         17,
	                                         3);
	   ServiceFlow UlServiceFlowUgs = wimax.CreateServiceFlow (ServiceFlow::SF_DIRECTION_UP,
	                                                           ServiceFlow::SF_TYPE_UGS,
	                                                           UlClassifierUgs);

	   UlServiceFlowUgs.SetIsEnabled(true);
	   UlServiceFlowUgs.SetMaximumLatency(1);
	   UlServiceFlowUgs.SetCsSpecification (ServiceFlow::IPV4);
	   UlServiceFlowUgs.SetServiceSchedulingType (ServiceFlow::SF_TYPE_UGS);
	   UlServiceFlowUgs.SetTrafficPriority (3);

	   DlServiceFlowUgs.SetIsEnabled(true);
	   DlServiceFlowUgs.SetMaximumLatency(2);

	   DlServiceFlowUgs.SetCsSpecification (ServiceFlow::IPV4);
	   DlServiceFlowUgs.SetServiceSchedulingType (ServiceFlow::SF_TYPE_UGS);
	   DlServiceFlowUgs.SetTrafficPriority (3);

	   //ss->AddServiceFlow (UlServiceFlowUgs);
	   ss->AddServiceFlow (DlServiceFlowUgs);

	   ss = 0;
	   bs = 0;

	   //checking MAC parameters

	   uint32_t RTGTTG_max=100;

	   uint32_t rtg=1;
	   uint32_t ttg=1;

	   float dis=m_distance;
	   float tProp=dis/300; //propagation time in us


	   //RTG = RTG_max;%min(2*Tprop, RTG_max);
	   //TTG = TTG_max;%min(2*Tprop, TTG_max);
	   if((2*tProp) > RTGTTG_max){
	   	 	rtg=RTGTTG_max;
	   	 	ttg=RTGTTG_max;
	   }else{
	   	 	rtg=2*tProp;
	   	 	ttg=2*tProp;
	   }

	     //N_RANGING_OPS per frame;
	     double ranging=1;

	     tProp=dis/300000000;
		 double fs= 144*10000000/125;
		 double tofdm=(1.25)/fs;
		 tofdm=tofdm*256;
		 ranging=(2*tProp)/tofdm;

		 if(ranging<3){
		  ranging=3;
		 }else{

		  ranging=floor(ranging);
		 }

	     Ptr<NetDevice> n;
	     Ptr<WimaxNetDevice> wd;
	     Ptr<BaseStationNetDevice> wbs;
	     n = bsDev.Get(0);
	     wbs=n->GetObject<BaseStationNetDevice>();
	     wd = n->GetObject<WimaxNetDevice>();

	   	 wd->SetRtg(rtg);
	   	 wd->SetTtg(ttg);
	   	 wbs->SetRangReqOppSize(ranging);

	   	std::cout << "\nMAC parameters: -------------------------------------------------------";
	   	std::cout << "\nRTG: " << wd->GetRtg() << std::endl;
	   	std::cout << "TTG: " << wd->GetTtg() << std::endl;
	   	std::cout << "Ranging: " << (uint32_t)wbs->GetRangReqOppSize() << std::endl;

	   NS_LOG_INFO ("Starting simulation.....");
	   Simulator::Run ();

	   Simulator::Destroy ();

	   std::cout << "\nSat. Throughput: " << ( (m_RxBytesDl *8)/m_runApplicationSeconds) << " bps"<<std::endl;
	   m_throughput=( (m_RxBytesDl *8)/m_runApplicationSeconds);

	   m_delay=1000*m_delayDl/m_RxPacketsDl;
	   std::cout << "\nDelay: " << m_delay << " ms " <<std::endl;

  }

