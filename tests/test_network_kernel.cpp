#include <fstream>
#include <regex>

#include <pops/network_kernel.hpp>

using namespace pops;

int test_node_status_at(
    const Network<int>& network, int row, int col, unsigned num_nodes_expected)
{
    int ret = 0;
    if (!network.has_node_at(row, col)) {
        std::cerr << "Node at " << row << ", " << col << " is missing\n";
        ret += 1;
    }
    auto nodes = network.get_nodes_at(row, col);
    if (nodes.size() != num_nodes_expected) {
        for (auto node : nodes) {
            std::cerr << "Node: " << node << "\n";
        }
        std::cerr << "There should be " << num_nodes_expected << " not " << nodes.size()
                  << " nodes at " << row << ", " << col << "\n";
        ret += 1;
    }
    return ret;
}

int test_create_network()
{
    int ret = 0;
    BBox<double> bbox;
    bbox.north = 40;
    bbox.south = 30;
    bbox.east = -70;
    bbox.west = -80;
    Network<int> network{
        bbox,
        0.01,
        0.01,
    };
    if (network.has_node_at(1, 1)) {
        std::cerr << "Empty network should not have a node at any cell\n";
        ret += 1;
    }
    std::stringstream node_stream{
        "1,-79.94,37.27\n"
        "2,-79.93,37.27\n"
        "3,-79.90,37.36\n"
        "4,-79.94,37.27\n"
        "5,-80.01,37.27\n"};
    std::stringstream segment_stream{
        "1,2,-79.93722002;37.27040997;-79.93691389;37.2708129;-79.93644276;37.27121265;-79.93600802;37.27151417;-79.9355576;37.27180514;-79.93486518;37.27225486;-79.93429095;37.272487\n"
        "3,4,-79.90224451;37.36736379;-79.90305243;37.36681812;-79.90349717;37.36651774;-79.90425106;37.36600856;-79.90500495;37.36549938;-79.90575884;37.3649902;-79.90647643;37.36450554;-79.90719401;37.36402089;-79.9079116;37.36353623;-79.90862919;37.36305157;-79.90934677;37.36256691;-79.91006436;37.36208225;-79.91078194;37.36159759;-79.91132198;37.361278;-79.91180896;37.3609898;-79.9125684;37.36058659;-79.91334929;37.36017879;-79.91411262;37.35978016;-79.91487594;37.35938153;-79.91544383;37.35904228;-79.91596933;37.35864311;-79.91642308;37.35822329;-79.91691486;37.3575987;-79.91745394;37.35677859;-79.91782659;37.3561825;-79.91834583;37.35535191;-79.91884291;37.35455677;-79.91933999;37.35376163;-79.91983707;37.35296648;-79.92031775;37.35219758;-79.92061594;37.35172059;-79.92093494;37.35114359;-79.92119894;37.35053759;-79.92139396;37.34992169;-79.92158898;37.34930578;-79.92174237;37.34882135;-79.92189576;37.34833691;-79.9221575;37.34751031;-79.92241923;37.34668371;-79.9226456;37.34596878;-79.92287198;37.34525384;-79.92309836;37.3445389;-79.92337425;37.34366757;-79.92365015;37.34279625;-79.92392605;37.34192492;-79.92420194;37.3410536;-79.92446494;37.34004289;-79.92462494;37.3392196;-79.92470894;37.3385566;-79.9248085;37.33763062;-79.92488055;37.33696046;-79.92494927;37.33632134;-79.92503299;37.33554262;-79.92514094;37.3345386;-79.92514194;37.3335916;-79.92506213;37.33255233;-79.92502013;37.33200544;-79.92496253;37.33125552;-79.92490494;37.3305056;-79.92494412;37.32978271;-79.92507694;37.3291406;-79.92533033;37.32830495;-79.92558372;37.3274693;-79.92580911;37.32672598;-79.9260345;37.32598266;-79.92625989;37.32523935;-79.92648528;37.32449603;-79.92683636;37.3233382;-79.92705436;37.32261925;-79.92725535;37.32201603;-79.92748889;37.32145503;-79.92801794;37.3203916;-79.92826894;37.3194926;-79.92837355;37.31886323;-79.92850671;37.31806213;-79.92863987;37.31726102;-79.92877302;37.31645992;-79.92890618;37.31565882;-79.92900694;37.3150526;-79.92912644;37.3145006;-79.92924594;37.3139486;-79.92942194;37.3132626;-79.92959794;37.3125766;-79.92970294;37.3118986;-79.92971221;37.31120799;-79.92965794;37.3103532;-79.92962941;37.309802;-79.92958951;37.3090309;-79.9295496;37.30825981;-79.92950969;37.30748871;-79.92946978;37.30671761;-79.92944281;37.30619641;-79.92941583;37.3056752;-79.92938448;37.30506931;-79.92935312;37.30446342;-79.92930706;37.30357348;-79.929261;37.30268354;-79.92921494;37.3017936;-79.92908461;37.30053767;-79.92898366;37.29979898;-79.92890638;37.29928835;-79.92866606;37.2981907;-79.92849844;37.29743755;-79.92833082;37.2966844;-79.9281632;37.29593126;-79.92801594;37.29526961;-79.92784894;37.29446261;-79.92768194;37.29365561;-79.92758617;37.2930236;-79.9274904;37.29239159;-79.92739742;37.29177799;-79.92729518;37.2911033;-79.92719294;37.29042861;-79.92707594;37.28945561;-79.92708294;37.28891861;-79.92711294;37.28835211;-79.92714294;37.28778561;-79.92719895;37.28685711;-79.92724929;37.28602273;-79.92729962;37.28518834;-79.92733072;37.28467284;-79.92736294;37.28413861;-79.92746594;37.28339061;-79.92761343;37.2825764;-79.92779244;37.28158823;-79.92788294;37.28108861;-79.92804276;37.28052758;-79.92820257;37.27996655;-79.92836807;37.27938552;-79.92853358;37.2788045;-79.92882894;37.27776761;-79.92898894;37.27729061;-79.92935094;37.27666361;-79.92991994;37.27610261;-79.93044294;37.27576661;-79.93101202;37.27547404;-79.9315811;37.27518147;-79.93207294;37.27492861;-79.93257014;37.27468982;-79.93323294;37.27441861;-79.93380594;37.27425961;-79.93429341;37.27414566;-79.93529894;37.27391061;-79.93577752;37.27374062;-79.93641199;37.27351525;-79.93708521;37.27333275;-79.93799896;37.27324423;-79.93892468;37.27325682;-79.93943233;37.27326595;-79.9399867;37.27327593;-79.94053063;37.27328965;-79.94134124;37.2733559\n"};
    network.load(node_stream, segment_stream);
    if (network.has_node_at(3727, -8001)) {
        std::cerr << "Node outside of bounding box was not ignored\n";
        ret += 1;
    }
    ret += test_node_status_at(network, 3727, -7994, 2);
    ret += test_node_status_at(network, 3727, -7993, 1);
    ret += test_node_status_at(network, 3736, -7990, 1);

    return ret;
}

using RawConfig = std::map<std::string, std::string>;

template<typename Stream>
RawConfig read_config(Stream& stream)
{
    RawConfig config;
    std::string line;
    while (std::getline(stream, line)) {
        std::regex delimeter(R"([\s]*:[\s]*)");
        std::smatch match;
        if (regex_search(line, match, delimeter)) {
            config[match.prefix()] = match.suffix();
        }
        else {
            throw std::runtime_error(std::string("Incorrect format at line: ") + line);
        }
    }
    return config;
}

BBox<double> bbox_from_config(const RawConfig& config)
{
    BBox<double> bbox;
    try {
        bbox.north = std::stod(config.at("north"));
        bbox.south = std::stod(config.at("south"));
        bbox.east = std::stod(config.at("east"));
        bbox.west = std::stod(config.at("west"));
    }
    catch (std::out_of_range&) {
        throw std::runtime_error(
            "At least one of north, south, east, west keys is missing");
    }
    return bbox;
}

int create_network_from_files(int argc, char** argv)
{
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0]
                  << " read|stats CONFIG_FILE NODE_FILE SEGMENT_FILE\n";
        return 1;
    }
    std::string command = argv[1];
    bool show_stats = false;
    if (command == "stats") {
        show_stats = true;
    }
    else if (command != "read") {
        std::cerr << "Unknown sub-command: " << command << "\n";
        std::cerr << "Supported sub-commands are: read, stats\n";
        return 1;
    }

    std::string config_file{argv[2]};
    std::ifstream config_stream{config_file};
    if (!config_stream.is_open()) {
        std::cerr << "Failed to open config file: " << config_file << "\n";
        return 1;
    }
    std::string node_file{argv[3]};
    std::ifstream node_stream{node_file};
    if (!node_stream.is_open()) {
        std::cerr << "Failed to open node file: " << node_file << "\n";
        return 1;
    }
    std::string segment_file{argv[4]};
    std::ifstream segment_stream{segment_file};
    if (!segment_stream.is_open()) {
        std::cerr << "Failed to open stream file: " << segment_file << "\n";
        return 1;
    }

    RawConfig config = read_config(config_stream);
    BBox<double> bbox = bbox_from_config(config);
    double nsres = std::stod(config["nsres"]);
    double ewres = std::stod(config["ewres"]);

    Network<int> network(bbox, nsres, ewres);
    network.load(node_stream, segment_stream);

    if (show_stats) {
        auto stats = network.collect_stats();
        for (const auto& item : stats) {
            std::cout << item.first << ": " << item.second << "\n";
        }
    }

    return 0;
}

int run_tests()
{
    int ret = 0;

    ret += test_create_network();

    std::cout << "Number of errors in the network kernel test: " << ret << "\n";
    return ret;
}

int main(int argc, char** argv)
{
    if (argc > 1)
        return create_network_from_files(argc, argv);
    else
        return run_tests();
}
