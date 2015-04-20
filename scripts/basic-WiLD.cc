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
This script simulate a WiFi Long Distance link (WiLD link) in a simple manner. It was developed for obtain throughput, latency, jitter and losses values when the link is working in the staturation point for different distances with a specific link configuration (MCS, Packet Size, Aggregation Threshold and AckTimeout Adjustment)


*      +-----+				+------------+		*
*      | AP  |      (<------------>)  	|    STA     |		*
*      +-----+          		+------------+		*
*      10.1.1.1          		   10.1.1.2		*


The accepted paramterers are:

	--mcs:       Modulation and codification scheme used. Different values for MCS can be used, 
		     from 0 to 15, which includes SISO and MIMO.
	--pSize:     Fixed packet size used for the simulation in bytes. This value correspond to APP 
		     level size. It is necessary to consider the protocol overhead since results are given at 
		     MAC level.
	--maxDist:   Maximum distance in meters to be tested. The script will test different distances 
  		     until reach this distance.
	--stepDist:  Step distance in meters to be tested. This parameter set the increase of distance
		     added in each simulation. For example, if "maxDist=30000" and "stepDist=5000", 
		     7 simultations will be performed at 5, 5005, 10005, 15005, 20005, 25005 and 30005 meters
	--aggTh:     Aggregation threshold used. This will set up the level of aggregation present in 
		     the WiLD link.
	--ackOpt:    Adjust ackTimeout for normal behavior. If true, this option will adjust the AckTimeout
		     optimally to get the maximum throughput in each distance. If false, the AckTimeout value
		     defined in the standard will be used, with the corresponding drop of performance starting
		     from a specific distance (usually 27 Km)[1]
	--pcap:      Generate Pcap trace file. If true, generate a Pcap file.

The following default values are used in this scripts:
	Standard:			IEEE 802.11n
	Frequency:			5 GHz
	Bandwidth: 		        20 MHz
	Guard Interval: 		Long (800ns)
	SlotTime:			Optimized according to [2]
	Fragmentation:			Disabled
	CTS/RTS:			Disabled
	BlockAck:			Disabled
	QoS:				EDCA not used, but available
	Transport Protocol:		UDP
	Propagation Model:		FixedRssLossModel and ConstantSpeedPropagationDelayModel
	Bidirectional flows:	Yes

Example of use:

./waf --run "scratch/basic-WiLD --mcs=0 --pSize=1300 --maxDist=30000 --stepDist=5000 --aggTh=0 --ackOpt=false --pcap=false"

In this experiment, results of throughput, latency, jitter and losses in the saturation point are given for different distances, from 0 Km to 30 Km in steps of 5 Km. The MCS used is 0, no aggregation is performed and the packet size at APP level is 1300 bytes. The ackOpt is no enabled, so a dramatic drop of the performance is expected at 27 Km. No pcap files are generated
*/

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/applications-module.h"
#include "ns3/wifi-module.h"
#include "ns3/mobility-module.h"
#include "ns3/internet-module.h"
#include "ns3/flow-monitor-module.h"
#include "ns3/gnuplot.h"


using namespace ns3;

NS_LOG_COMPONENT_DEFINE ("TUCAN3G");


class Experiment
{
	public:
		Experiment ();
		void RunSim (uint32_t dataRate, uint32_t distance, uint32_t packetSize, std::string phyMode,
				uint32_t aggregationThreshold, float acktimeout,uint32_t slot_time, uint32_t runApplicationSeconds, int32_t rssForFixedRssLossModel,
				 bool pcapTrace);
		uint64_t GetThroughput ();
		uint64_t GetDelay ();
		uint64_t GetJitter ();
		uint64_t GetLosses ();

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
		std::string m_phyMode;
		uint32_t m_aggregationThreshold;
		uint32_t m_throughput;
		uint32_t m_delay;
		uint32_t m_jitter;
		uint32_t m_losses;
		uint32_t m_runApplicationSeconds;
		int32_t m_rssi;
		float m_acktimeout;
		uint32_t m_slot_time;
};


Experiment::Experiment () :
		port1(4002),
		port2(4000),

		m_throughput(0)
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
	 std::cout <<"   --mcs:       Modulation and codification scheme used [0]\n";
	 std::cout <<"   --pSize:     Fixed packet size used for the simulation in bytes [1372]\n";
	 std::cout <<"   --maxDist:   Maximum distance in meters to be tested [50000]\n";
	 std::cout <<"   --stepDist:  Step distance in meters to be tested [5000]\n";
	 std::cout <<"   --aggTh:     Aggregation threshold used [8192]\n";
	 std::cout <<"   --olRes:     Offered Load resolution [10]\n";
         std::cout <<"   --pcap:      Generate Pcap trace file (boolean) [false]\n";
}

  int main (int argc, char *argv[])
  {

  	Experiment experiment;

  	// Command line arguments
  	int32_t mcs = 0;	//mcs 0-15
  	int32_t aggregationThreshold =8192;	//aggregation threshold
  	int32_t maxDistance=50000; //max distance in meters
  	int32_t stepDistance=5000; //step distance in meters ( < maxDistance)
  	bool ackTimeoutOpt=true;	//adjust ackTimeout for optimal behavior
  	int32_t packetSize = 1372;	//Fixed packet size used for the simulation (app level)
        bool pcapTrace = false;
  	CommandLine cmd;
  	cmd.AddValue ("mcs", "Modulation and codification scheme used", mcs);
  	cmd.AddValue ("pSize", "Fixed packet size used for the simulation in bytes", packetSize);
  	cmd.AddValue ("maxDist", "Maximum distance in meters to be tested",maxDistance);
  	cmd.AddValue ("stepDist", "Step distance in meters to be tested", stepDistance);
    cmd.AddValue ("aggTh", "Aggregation threshold used ( 0 < aggTh < 8192 )", aggregationThreshold);
  	cmd.AddValue ("ackOpt", "Adjust ackTimeout for normal behavior ", ackTimeoutOpt);
    cmd.AddValue ("pcap", "Create Pcap trace file (boolean)", pcapTrace);
  	cmd.Parse (argc, argv);
        
    std::cout << "\n******************************************************************";
  	std::cout << "\nBeging simulation with the following configuration:";

  	std::cout << "\nBandwidth -> 20MHz";
  	if(mcs==0){
  		std::cout << "\nMCS -> " << mcs << " (Using Default)";
  	}else{
  		std::cout << "\nMCS -> " << mcs;
  	}
  	if(aggregationThreshold==8192){
  	  	std::cout << "\nAggregation -> " << aggregationThreshold << " (Using Default)";
  	}else{
  	  	std::cout << "\nAggregation -> " << aggregationThreshold;
  	}
        if(packetSize==1372){
    	        std::cout << "\nPacketSize -> " << packetSize << " (Using Default)";
        }else{
    	        std::cout << "\nPacketSize -> " << packetSize;
        }
        if(maxDistance==50000){
                std::cout << "\nMaxDistance -> " << maxDistance << " (Using Default)";
        }else{
                 std::cout << "\nMaxDistance -> " << maxDistance;
        }
        if(stepDistance==5000){
                 std::cout << "\nStepDistance -> " << stepDistance << " (Using Default)";
        }else{
                 std::cout << "\nStepDistance -> " << stepDistance;
        }
        std::cout << "\n******************************************************************\n";

  	uint32_t slot_init=9; //ns
  	uint32_t slot=slot_init;

  	uint32_t distance_init = 5;
  	uint32_t distance=distance_init;

  	int32_t rssi = 0;
  	uint32_t j = 0;
  	float acktimeout = 0.0002; //will be overwritten

  	if((mcs < 0) || (mcs > 15)){
  		NS_LOG_UNCOND ("MCS has to be one of these: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 ");
  		PrintHelp();
  		return 0;
  	}else{
  		  j = mcs;
  	}

  	if(packetSize <= 0){
  	  	NS_LOG_UNCOND ("The packet size must be positive");
  	    PrintHelp();
  	  	return 0;
  	}if(maxDistance < 0){
  	  	NS_LOG_UNCOND ("The maximum distance must be positive");
  	    PrintHelp();
  	  	return 0;
  	}if((stepDistance < 0) ||(stepDistance > maxDistance)) {
  	  	NS_LOG_UNCOND ("The step distance must be positive and lower than the maximum distance");
  	    PrintHelp();
  	  	return 0;
  	}if((aggregationThreshold < 0) ||(aggregationThreshold > 8192)) {
  	  	NS_LOG_UNCOND ("The maximum aggregation threshold allowed is 8192");
  	    PrintHelp();
  	  	return 0;
  	}

  	 //simulation parameters

  	 uint32_t runApplicationSeconds = 15;

  	//physical data rates MCS0 - MCS15
  	std::string phyModeL[16] = {"OfdmRate6_5MbpsBW20MHz",
                   "OfdmRate13MbpsBW20MHz",
                   "OfdmRate19_5MbpsBW20MHz",
                   "OfdmRate26MbpsBW20MHz",
                   "OfdmRate39MbpsBW20MHz",
                   "OfdmRate52MbpsBW20MHz",
                   "OfdmRate58_5MbpsBW20MHz",
                   "OfdmRate65MbpsBW20MHz",
                   "OfdmRate13MbpsBW20MHz",
                   "OfdmRate26MbpsBW20MHz",
                   "OfdmRate39MbpsBW20MHz",
                   "OfdmRate52MbpsBW20MHz",
                   "OfdmRate1_5MbpsBW5MHz", 	        //equivalent to 78Mbps
                   "OfdmRate2_25MbpsBW5MHz", 	        //equivalent to 104Mbps
                   "OfdmRate3MbpsBW5MHz",		//equivalent to 117Mbps
                   "OfdmRate4_5MbpsBW5MHz"};		//equivalent to 130Mbps

  	// injected datarates in each node
  	uint32_t dataRateL[16]={3250,6500,9750,13000,19500,26000,29250,32500,
  							6500,13000,19500,26000,39000,52000,58500,65000};

        
    std::string phyModeL_printable[16]= phyModeL;
  	phyModeL_printable[12]="OfdmRate78MbpsBW20MHz";
  	phyModeL_printable[13]="OfdmRate104MbpsBW20MHz";
  	phyModeL_printable[14]="OfdmRate114MbpsBW20MHz";
  	phyModeL_printable[15]="OfdmRate130MbpsBW20MHz";

    std::stringstream aa;
	aa << aggregationThreshold << "bytes";
	std::string thres = aa.str();
       
  	std::string file_delay  = "802.11n_" + phyModeL_printable[j] + "Agregagc"+ thres +"_delay.ods";
	std::string file_throughput  = "802.11n_" + phyModeL_printable[j] + "Agregagc"+ thres +"_throguhput.ods";
	std::string file_jitter  = "802.11n_" + phyModeL_printable[j] + "Agregagc"+ thres +"_jitter.ods";
	std::string file_losses  = "802.11n_" + phyModeL_printable[j] + "Agregagc"+ thres +"_losses.ods";


    //preparing output files
  	std::ofstream out_th(file_throughput.c_str(), std::ios::app);
  	out_th << "Distance (m) , Sat. Throughput (bps)"<<	std::ends;
  	out_th << std::endl;
  	out_th.close();
  	std::ofstream out_delay(file_delay.c_str(), std::ios::app);
  	out_delay << "Distance (m) , Sat. Delay (us)"<<	std::ends;
  	out_delay << std::endl;
  	out_delay.close();
  	std::ofstream out_jitter(file_jitter.c_str(), std::ios::app);
  	out_jitter  << "Distance (m) , Sat. Jitter (us)"<<	std::ends;
	out_jitter << std::endl;
	out_jitter.close();
	std::ofstream out_losses(file_losses.c_str(), std::ios::app);
	out_losses << "Distance (m) , Sat. Losses "<<	std::ends;
  	out_losses << std::endl;
  	out_losses.close();

  	slot=(slot_init-2)+(2*(distance/1000)/0.3);
  	if(slot < slot_init){
  	    slot=slot_init;
  	}
  	double tprop=0;
  	if(ackTimeoutOpt){
  		tprop=distance / 300.0;
  	}else{
  		tprop=93; //if coverage class 31 IEEE 802.11-2012 Table 8-56. MAX 27 Km

  	}
  	acktimeout=(10+52+2*tprop+20);
  	while(distance < (maxDistance+1)){
  			std::cout << "\nSimulating with: -------------------------------------------------------";
  			std::cout << "\nSlotTime -> " << slot;
  			std::cout << "\nAggregation -> " << aggregationThreshold;
  			std::cout << "\nPacketSize -> " << packetSize;
  			std::cout << "\nDistance -> " << distance;
  			std::cout << "\nMaxDistance -> " << maxDistance;
  			std::cout << "\nInjected datarate -> " << dataRateL[j]*2;

  		    experiment = Experiment();
  		    experiment.RunSim(dataRateL[j], distance, packetSize, phyModeL[j],
  				    aggregationThreshold, acktimeout,slot,
  				    runApplicationSeconds,rssi, pcapTrace);

  	       //fill output files with output values
               std::stringstream ss;
               ss <<distance;
               std::string sdistance = ss.str();
  	       std::ofstream out_th2(file_throughput.c_str(), std::ios::app);
  	       out_th2 << sdistance <<" ,"<< experiment.GetThroughput() << " ,"<<	std::ends;
  	       out_th2 << std::endl;
  	       out_th2.close();
  	       std::ofstream out_delay2(file_delay.c_str(), std::ios::app);
  	       out_delay2 << sdistance <<" ,"<< experiment.GetDelay() << " ,"<<	std::ends;
  	       out_delay2 << std::endl;
  	       out_delay2.close();
  	       std::ofstream out_jitter2(file_jitter.c_str(), std::ios::app);
  	       out_jitter2 << sdistance <<" ,"<< experiment.GetJitter() << " ,"<<	std::ends;
  	       out_jitter2 << std::endl;
  	       out_jitter2.close();
  	       std::ofstream out_losses2(file_losses.c_str(), std::ios::app);
  	       out_losses2 << sdistance <<" ,"<< experiment.GetLosses() << " ,"<<	std::ends;
  	       out_losses2 << std::endl;
  	       out_losses2.close();

  	       distance=distance+stepDistance;
  	       slot=(slot_init-2)+(2*(distance/1000)/0.3);
  	       if(slot < slot_init){
  	    	      slot=slot_init;
  	       }
  	       if(ackTimeoutOpt){
                   tprop=distance / 300.0;
           }else{
                   tprop=93; //if coverage class 31 IEEE 802.11-2012 Table 8-56 . MAX 27 Km
           }
           acktimeout=(10+52+2*tprop+20);
  	}
  	NS_LOG_UNCOND ("Done... ");
        std::cout << "\n ---------------------------------------------";
          	std::cout << "\nThroughput values are in: \n" << file_throughput;
          	std::cout << "\nDelay values are in: \n" << file_delay;
          	std::cout << "\nJitter values are in: \n" << file_jitter;
          	std::cout << "\nLosses values are in: \n" << file_losses;
        if(pcapTrace){
          	std::cout << "\nDL Pcap traces are in: wifi-simple-infra-0-0.pcap";
          	std::cout << "\nUL Pcap traces are in: wifi-simple-infra-1-0.pcap \n";
        }
  	return 0;
  }

  uint64_t Experiment::GetThroughput (){
	 return m_throughput;
  }
  uint64_t Experiment::GetDelay (){
  	 return m_delay;
    }
  uint64_t Experiment::GetJitter (){
  	 return m_jitter;
    }
  uint64_t Experiment::GetLosses (){
  	 return m_losses;
    }

 void
 Experiment::RunSim (uint32_t dataRate, uint32_t distance, uint32_t packetSize, std::string phyMode,
		  uint32_t aggregationThreshold, float acktimeout, uint32_t slot, uint32_t runApplicationSeconds,
			int32_t rssi, bool pcapTrace)
  {

	  m_dataRate = dataRate;
	  m_distance = distance;
	  m_packetSize = packetSize;
	  m_phyMode = phyMode;
	  m_aggregationThreshold = aggregationThreshold;

	  m_acktimeout = acktimeout;

	  m_runApplicationSeconds = runApplicationSeconds;

	  m_rssi = rssi;
	  m_slot_time = slot;

	  // ad-hoc mode
  	  NodeContainer tas;
  	  tas.Create (2);

  	  YansWifiChannelHelper channel = YansWifiChannelHelper::Default ();
  	  YansWifiPhyHelper phy = YansWifiPhyHelper::Default ();
  	  phy.SetPcapDataLinkType (YansWifiPhyHelper::DLT_IEEE802_11_RADIO);

  	  channel.SetPropagationDelay ("ns3::ConstantSpeedPropagationDelayModel", "Speed", DoubleValue (300000000.0) );
  	  channel.AddPropagationLoss ("ns3::FixedRssLossModel", "Rss", DoubleValue (m_rssi)); //fixed received signal


  	  phy.SetChannel (channel.Create ());
  	  phy.Set ("ShortGuardEnabled",BooleanValue(false));

  	  // turn off RTS/CTS for frames below 15200 bytes. This mean practically that always is turn off
  	  Config::SetDefault ("ns3::WifiRemoteStationManager::RtsCtsThreshold", UintegerValue (15200));
  	  // disable fragmentation for frames below 15200 bytes. This mean practically that always is disabled
  	  Config::SetDefault ("ns3::WifiRemoteStationManager::FragmentationThreshold", UintegerValue (15200));
  	  // Fix non-unicast data rate to be the same as that of unicast
  	  Config::SetDefault ("ns3::WifiRemoteStationManager::NonUnicastMode",StringValue (phyMode));

  	  WifiHelper wifi = WifiHelper::Default ();

  	  wifi.SetStandard (WIFI_PHY_STANDARD_80211n_2_4GHZ);
  	  //this method is not relevant since the main mac parameters such us DIFS,RIFS, TimeSlot,...
  	  //are preconfigured manually

  	  wifi.SetRemoteStationManager ("ns3::ConstantRateWifiManager",
  	                                  "DataMode",StringValue (phyMode),
  	                                  "ControlMode",StringValue (phyMode));

      // ad-hoc mode
  	  Ssid ssid = Ssid ("ns-3-ssid");
  	  HtWifiMacHelper wifiMac = HtWifiMacHelper::Default ();

  	  wifiMac.SetType ("ns3::AdhocWifiMac",
  	    		"Ssid", SsidValue (ssid),
  	    								"Sifs", TimeValue (MicroSeconds (10)),
  	    								"Slot", TimeValue (MicroSeconds (m_slot_time)),
  	    								"EifsNoDifs", TimeValue (MicroSeconds (10+304)),
  	    								"Pifs", TimeValue (MicroSeconds (10+m_slot_time)),
  	    		      	                "AckTimeout", TimeValue (MicroSeconds (m_acktimeout)),
  	    		      	                "CtsTimeout", TimeValue (MicroSeconds (m_acktimeout)));
  	    		      	                //"MaxPropagationDelay", TimeValue (Seconds (100000.0 / 300000000.0)));

  	NetDeviceContainer devices = wifi.Install (phy, wifiMac, tas);


  	wifiMac.SetMsduAggregatorForAc (AC_BE, "ns3::MsduStandardAggregator", "MaxAmsduSize", UintegerValue (aggregationThreshold));
  	wifiMac.SetBlockAckThresholdForAc (AC_BE, 0);


  	 //Mobility
  	 // Note that with FixedRssLossModel, the positions below are not
  	 // used for received signal strength.
  	  MobilityHelper mobility;
  	  mobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  	  Ptr<ListPositionAllocator> positionAlloc = CreateObject<ListPositionAllocator> ();
  	  positionAlloc->Add (Vector (0.0, 0.0, 0.0));
  	  positionAlloc->Add (Vector (distance, 0.0, 0.0));
  	  mobility.SetPositionAllocator (positionAlloc);

  	  mobility.SetMobilityModel ("ns3::ConstantPositionMobilityModel");
  	  mobility.Install (tas);

  	  InternetStackHelper stack;
  	  stack.Install (tas);

  	  Ipv4AddressHelper ipv4; //Assign IP Addresses
  	  ipv4.SetBase ("10.1.1.0", "255.255.255.0");
  	  Ipv4InterfaceContainer i = ipv4.Assign (devices);

  	  GlobalValue::Bind ("ChecksumEnabled", BooleanValue (true));

  	  Ipv4GlobalRoutingHelper::PopulateRoutingTables ();
  	  PopulateArpCache ();

  	  std::stringstream ss1;
  	  ss1 << m_dataRate;
  	  std::string rate = ss1.str();
  	  rate = rate + "kbps";

  	  //Flow 0 (10.1.1.2 -> 10.1.1.1)

  	  OnOffHelper onoff1 ("ns3::UdpSocketFactory", InetSocketAddress (i.GetAddress (1), port1));
  	  onoff1.SetAttribute ("OnTime", StringValue ("ns3::ConstantRandomVariable[Constant=1]"));
  	  onoff1.SetAttribute ("OffTime", StringValue ("ns3::ConstantRandomVariable[Constant=0]"));
  	  onoff1.SetAttribute ("DataRate", DataRateValue (DataRate (rate)));
  	  onoff1.SetAttribute ("PacketSize", UintegerValue (m_packetSize));
  	  //onoff1.SetAttribute("AccessClass",UintegerValue (1));
  	  AddressValue remoteAddress1 (InetSocketAddress (i.GetAddress (0), port1));
  	  onoff1.SetAttribute ("Remote", remoteAddress1);

  	  ApplicationContainer temp1 = onoff1.Install (tas.Get(1));

  	  temp1.Start(Seconds (1));//avoid collisions at the begining 
  	  temp1.Stop (Seconds (m_runApplicationSeconds+1));

  	  Ptr<Socket> sink1 = SetupPacketReceive1 (i.GetAddress (0), tas.Get(0));

  	  //Flow 1 (10.1.1.1 -> 10.1.1.2)

  	  OnOffHelper onoff2 ("ns3::UdpSocketFactory", InetSocketAddress (i.GetAddress (0), port2));
  	  onoff2.SetAttribute ("OnTime", StringValue ("ns3::ConstantRandomVariable[Constant=1]"));
  	  onoff2.SetAttribute ("OffTime", StringValue ("ns3::ConstantRandomVariable[Constant=0]"));
  	  onoff2.SetAttribute ("DataRate", DataRateValue (DataRate (rate)));
  	  onoff2.SetAttribute ("PacketSize", UintegerValue (m_packetSize));
  	  //onoff2.SetAttribute("AccessClass",UintegerValue (1));
  	  AddressValue remoteAddress2 (InetSocketAddress (i.GetAddress (1), port2));
  	  onoff2.SetAttribute ("Remote", remoteAddress2);

  	  ApplicationContainer temp2 = onoff2.Install (tas.Get(0));
  	  temp2.Start(Seconds (1.00001));
  	  temp2.Stop (Seconds (m_runApplicationSeconds+1));

  	  Ptr<Socket> sink2 = SetupPacketReceive2 (i.GetAddress (1), tas.Get(1));


  	  if(pcapTrace){
  		  phy.EnablePcap ("wifi-simple-infra", devices);
  	  }


  	  FlowMonitorHelper flowmon;
  	  Ptr<FlowMonitor> monitor = flowmon.InstallAll ();



  	  Ptr<NetDevice> n1,n2;
  	  Ptr<WifiNetDevice> wd1,wd2;
  	  n1 = devices.Get(0);
  	  wd1 = n1->GetObject<WifiNetDevice>();   // This was the secret sauce to getting access to that mac object
  	  n2 = devices.Get(1);
  	  wd2 = n2->GetObject<WifiNetDevice>();

  	  PointerValue ptr1;
  	  wd1->GetMac()->GetAttribute("DcaTxop", ptr1);
  	  Ptr<DcaTxop> dca1 = ptr1.Get<DcaTxop>();
  	  //dca1->SetAifsn(3);

  	  std::cout << "\n\nChecking MAC Parameters:" << std::endl;
  	  std::cout << "SlotTime: " << wd1->GetMac()->GetSlot() << std::endl; // Prints out default
  	  std::cout << "AckTimeout: " << wd1->GetMac()->GetAckTimeout() << std::endl;
  	  std::cout << "SIFS: " << wd1->GetMac()->GetSifs() << std::endl;
  	  std::cout << "RIFS: " << wd1->GetMac()->GetRifs() << std::endl;
  	  std::cout << "EIFS " << wd1->GetMac()->GetEifsNoDifs() << std::endl;
  	  std::cout << "AIFSN " << dca1->GetAifsn() << std::endl;
  	  std::cout << "PIFS " << wd1->GetMac()->GetPifs() <<"\n"<<std::endl;

  	  Simulator::Stop (Seconds (1+m_runApplicationSeconds+1));
  	  Simulator::Run ();

  	  uint32_t throughput1=0;
  	  uint32_t throughput2=0;
  	  uint64_t delay1=0;
  	  uint64_t delay2=0;
  	  uint64_t jitter1=0;
  	  uint64_t jitter2=0;
  	  uint32_t losses1=0;
  	  uint32_t losses2=0;
  	  uint32_t a=0;
  	  uint32_t b=0;

  	  //configuring flowmonitor
  	  double Throughput =0;
  	  monitor->CheckForLostPackets ();
  	  		    Ptr<Ipv4FlowClassifier> classifier = DynamicCast<Ipv4FlowClassifier> (flowmon.GetClassifier ());
  	  		        std::map<FlowId, FlowMonitor::FlowStats> stats = monitor->GetFlowStats ();
  	  		    for (std::map<FlowId, FlowMonitor::FlowStats>::const_iterator i = stats.begin (); i != stats.end (); ++i)
  	  		    {
  	  		          Ipv4FlowClassifier::FiveTuple t = classifier->FindFlow (i->first);
  	  		          std::cout << "Flow " << i->first - 1 << " (" << t.sourceAddress << " -> " << t.destinationAddress << ")\n";
  	  		          //Throughput = i->second.rxBytes * 8.0 / m_runApplicationSeconds;
  	  		          Throughput = i->second.rxBytes * 8.0 / (i->second.timeLastRxPacket.GetSeconds()- i->second.timeFirstTxPacket.GetSeconds());

  	  		          std::cout << "Tiempo " << i->second.timeLastRxPacket.GetSeconds()- i->second.timeFirstTxPacket.GetSeconds() << " bps\n";
  	  		          std::cout << "Tiempo2 " << m_runApplicationSeconds << " \n";
  	  		          if((i->first - 1) ==0){//flow 0

  	  		        	throughput1=Throughput*(m_packetSize+8+20+60)/(m_packetSize+20+8);
  	  		        	a=(i->second.rxPackets * (m_packetSize) * 8)/( m_runApplicationSeconds);
  	  		        	//std::cout << i->first - 1 << "  Tx Packets:   " << i->second.txPackets << "\n";
  	  		        	//std::cout << i->first - 1 << "  Rx Packets:   " << i->second.rxPackets << "\n";
  	  		        	//std::cout << i->first - 1 << "  Rx bytes:   " << i->second.rxBytes << "\n";
  	  		        	std::cout << "Throughput APP: " << a  << " bps\n";
  	  		        	std::cout << "Throughput MAC: " << throughput1 << " bps\n";
  	  		        	delay1 = i->second.delaySum.GetInteger() / i->second.rxPackets / 1000;
  	  		        	losses1=i->second.lostPackets;
  	  		        	jitter1=i->second.jitterSum.GetInteger() / i->second.rxPackets /1000;


  	  		          }else{ //flow1

  	  		        	b=(i->second.rxPackets * (m_packetSize) * 8)/( m_runApplicationSeconds);
  	  		        	throughput2=Throughput*(m_packetSize+8+20+60)/(m_packetSize+20+8);
  	  		        	//std::cout << i->first - 1 << "  Tx Packets:   " << i->second.txPackets << "\n";
  	  		        	//std::cout << i->first - 1 << "  Rx Packets:   " << i->second.rxPackets << "\n";
  	  		        	//std::cout << i->first - 1 << "  Rx bytes:   " << i->second.rxBytes << "\n";
  	  		        	std::cout << "Throughput APP: " << b  << " bps\n";
  	  		        	std::cout << "Throughput MAC/IP: " << throughput2 << " bps\n";
  	  		        	delay2 = i->second.delaySum.GetInteger() / i->second.rxPackets / 1000;
  	  		        	losses2=i->second.lostPackets;
  	  		        	jitter2=i->second.jitterSum.GetInteger() / i->second.rxPackets /1000;
  	  		          }

  	  		    }

  	  		//update values
  	  		m_delay=(delay1+delay2)*0.5;
  	  		m_jitter=(jitter1+jitter2)*0.5;
  	  		m_throughput=throughput1+throughput2;
  	  		m_losses=(losses1+losses2);

  	  	   std::cout << "****************************************************************************";
  	  		std::cout << "\nLatency: " << m_delay << " us\n";
  	  		std::cout << "Jitter: " << m_jitter << " us\n";
  	  		std::cout << "Losses: " << m_losses << " packets\n";
  	  		std::cout << "Total Throughput Mac: " << m_throughput << " bps\n";
  	  		std::cout << "Total Throughput App: " << a+b << " bps\n";
  	  	   std::cout << "****************************************************************************";

  	  Simulator::Destroy ();


  }

