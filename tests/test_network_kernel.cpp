#include <fstream>
#include <regex>
#include <random>

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

int compare_network_statistics(
    Network<int>::Statistics expected, Network<int>::Statistics actual)
{
    int ret{0};
    for (const auto& item : expected) {
        auto it = actual.find(item.first);
        if (it == actual.end()) {
            std::cerr << "Key " << item.first << " is missing in statistics\n";
            ++ret;
            continue;
        }
        auto value = it->second;
        if (value != item.second) {
            std::cerr << "Value for key " << item.first << " differs, expected "
                      << item.second << " got " << value << "\n";
            ++ret;
        }
    }
    return ret;
}

int test_travel_network()
{
    int ret = 0;
    BBox<double> bbox;
    bbox.north = 10;
    bbox.south = 0;
    bbox.east = 30;
    bbox.west = 20;
    Network<int> network{
        bbox,
        1,
        1,
        1,
    };
    std::stringstream node_stream{
        "1,22.0,7.0\n"
        "2,27.0,6.0\n"
        "3,21.0,1.0\n"
        "8,28.0,3.0\n"};
    std::stringstream segment_stream{
        "1,2,22.0;7.0;22.7;7.5;23.5;8.0;24.3;8.5;25.0;9.0;25.8;9.0;26.6;9.0;27.4;9.0;28.2;9.0;29.0;9.0;29.0;8.0;29.0;7.0;29.0;6.0;29.0;5.0;29.0;4.0;29.0;3.0;29.0;2.0;29.0;1.0;28.3;1.3;27.5;1.5;26.8;1.8;26.0;2.0;26.2;2.8;26.4;3.6;26.6;4.4;26.8;5.2;27.0;6.0\n"
        "3,8,21.0;1.0;21.2;1.9;21.3;2.7;21.4;3.6;21.6;4.4;21.7;5.3;21.9;6.1;22.0;7.0;22.8;6.8;23.7;6.7;24.5;6.5;25.3;6.3;26.2;6.1;27.0;6.0;27.3;5.2;27.5;4.5;27.8;3.7;28.0;3.0\n"};
    network.load(node_stream, segment_stream);

    std::random_device generator;
    double time = 9;
    int start_row = 9;
    int start_col = 1;
    if (!network.has_node_at(start_row, start_col)) {
        std::cerr << "Expected node at " << start_row << ", " << start_col << "\n";
        network.dump_yaml(std::cerr);
        ret += 1;
    }
    int end_row;
    int end_col;
    std::tie(end_row, end_col) = network.travel(start_row, start_col, time, generator);
    std::cerr << "from (" << start_row << ", " << start_col << ") to (" << end_row
              << ", " << end_col << ") in " << time << "\n";
    network.dump_yaml(std::cerr);
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
    Network<int> network{bbox, 0.01, 0.01, 42};
    if (network.has_node_at(1, 1)) {
        std::cerr << "Empty network should not have a node at any cell\n";
        ret += 1;
    }
    std::stringstream node_stream{
        "1,-79.937,37.270\n"
        "2,-79.934,37.272\n"
        "3,-79.902,37.367\n"
        "4,-79.941,37.273\n"
        "5,-80.015,37.279\n"};
    std::stringstream segment_stream{
        "1,2,-79.937;37.270;-79.936;37.270;-79.936;37.271;-79.936;37.271;-79.936;37.271;-79.934;37.272;-79.934;37.272\n"
        "3,4,-79.902;37.367;-79.903;37.366;-79.903;37.366;-79.904;37.366;-79.905;37.365;-79.905;37.36;-79.920;37.352;-79.93;37.273;-79.940;37.273;-79.941;37.273\n"};
    network.load(node_stream, segment_stream);
    int out_row;
    int out_col;
    std::tie(out_row, out_col) = network.xy_to_row_col(-80.015, 37.279);
    if (network.has_node_at(out_row, out_col)) {
        std::cerr << "Node outside of bounding box was not ignored\n";
        ret += 1;
    }
    ret += test_node_status_at(network, 272, 6, 2);
    ret += test_node_status_at(network, 263, 9, 1);
    ret += test_node_status_at(network, 272, 5, 1);

    Network<int>::Statistics expected_stats{
        {"num_nodes", 4},
        {"num_segments", 2},
        {"num_nodes_with_segments", 4},
        {"num_standalone_nodes", 0}};
    ret += compare_network_statistics(expected_stats, network.collect_statistics());

    if (ret) {
        network.dump_yaml(std::cerr);
    }

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
                  << " read|stats|write CONFIG_FILE NODE_FILE SEGMENT_FILE\n";
        return 1;
    }
    std::string command = argv[1];
    bool show_stats = false;
    bool write_network = false;
    if (command == "stats") {
        show_stats = true;
    }
    else if (command == "write") {
        write_network = true;
    }
    else if (command != "read") {
        std::cerr << "Unknown sub-command: " << command << "\n";
        std::cerr << "Supported sub-commands are: read, stats, write\n";
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

    Network<int> network(bbox, nsres, ewres, 42);
    network.load(node_stream, segment_stream);

    if (show_stats) {
        auto stats = network.collect_statistics();
        for (const auto& item : stats) {
            std::cout << item.first << ": " << item.second << "\n";
        }
    }
    if (write_network) {
        network.dump_yaml(std::cout);
    }

    return 0;
}

int run_tests()
{
    int ret = 0;

    ret += test_create_network();
    ret += test_travel_network();

    if (ret)
        std::cerr << "Number of errors in the network kernel test: " << ret << "\n";
    return ret;
}

int main(int argc, char** argv)
{
    if (argc > 1)
        return create_network_from_files(argc, argv);
    else
        return run_tests();
}
