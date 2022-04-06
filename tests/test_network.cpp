#include <fstream>
#include <regex>
#include <random>

#include <pops/network.hpp>

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

int test_bbox_functions()
{
    int ret = 0;
    BBox<double> bbox;
    bbox.north = 10;
    bbox.south = 0;
    bbox.east = 30;
    bbox.west = 20;
    Network<int> network{bbox, 1, 1};
    if (network.xy_out_of_bbox(25, 5)) {
        std::cerr << "XY 25, 5 should be in, not out\n";
        ++ret;
    }
    if (!network.xy_out_of_bbox(25, 11)) {
        std::cerr << "XY 25, 11 should be out, not in\n";
        ++ret;
    }
    if (network.row_col_out_of_bbox(2, 3)) {
        std::cerr << "Row 2, col 3 should be in, not out\n";
        ++ret;
    }
    if (!network.row_col_out_of_bbox(20, 3)) {
        std::cerr << "Row 20, col 3 should be out, not in\n";
        ++ret;
    }
    if (network.cell_out_of_bbox(std::make_pair(9, 8))) {
        std::cerr << "Cell 9, 8 should be in, not out\n";
        ++ret;
    }
    if (!network.cell_out_of_bbox(std::make_pair(-1, 3))) {
        std::cerr << "Cell -1, 3 should be out, not in\n";
        ++ret;
    }
    if (network.xy_out_of_bbox(bbox.east, bbox.north)) {
        std::cerr << "Bbox EN corner XY should be in, not out\n";
        ++ret;
    }
    if (network.xy_out_of_bbox(bbox.west, bbox.south)) {
        std::cerr << "Bbox WS corner XY should be in, not out\n";
        ++ret;
    }
    if (network.cell_out_of_bbox(network.xy_to_row_col(bbox.east, bbox.north))) {
        std::cerr << "Bbox EN corner cell should be in, not out\n";
        ++ret;
    }
    if (network.cell_out_of_bbox(network.xy_to_row_col(bbox.west, bbox.south))) {
        std::cerr << "Bbox WS corner cell should be in, not out\n";
        ++ret;
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
    Network<int> network{bbox, 1, 1};
    std::stringstream network_stream{
        "1,2,21.4;7.5;22.3;7.2\n"
        "1,4,21.4;7.5;21.9;8.0;22.5;8.6\n"
        "5,1,27.5;1.5;26.7;1.4;25.9;1.3;25.2;1.1;24.5;1.7;23.9;2.3;23.2;2.8;22.5;2.3;21.8;1.8;21.2;1.3;20.5;1.4;20.7;2.3;20.8;3.2;20.9;4.0;21.0;4.9;21.1;5.7;21.3;6.6;21.4;7.5\n"
        "2,8,22.3;7.2;23.2;7.1;24.0;6.9;24.8;6.7;25.7;6.6;26.5;6.4\n"
        "8,10,26.5;6.4;26.8;5.7;27.2;4.9;27.5;4.2;27.9;3.5;28.2;2.7\n"
        "11,8,28.3;9.0;28.5;8.1;28.6;7.3;28.2;6.8;27.7;6.3;27.1;6.4;26.5;6.4\n"
        "2,5,22.3;7.2;23.0;7.7;23.7;8.1;24.3;8.5;25.0;9.0;25.8;9.0;26.6;9.0;27.4;9.0;28.2;9.0;29.0;9.0;29.0;8.0;29.0;7.0;29.0;6.0;29.0;5.0;29.0;4.0;29.0;3.0;29.0;2.0;29.0;1.0;28.2;1.3;27.5;1.5\n"
        "5,8,27.5;1.5;26.7;1.8;26.0;2.0;26.1;2.9;26.2;3.8;26.3;4.7;26.4;5.5;26.5;6.4\n"};
    network.load(network_stream);

    std::default_random_engine generator;
    const int num_distances = 12;
    int current_distance = 0;
    std::array<int, num_distances> distances;
    std::generate_n(distances.begin(), num_distances, [&current_distance] {
        return ++current_distance;
    });
    std::array<std::pair<int, int>, num_distances> correct_destinations{
        {{8, 8},
         {9, 9},
         {8, 9},
         {7, 9},
         {6, 9},
         {5, 9},
         {4, 9},
         {3, 9},
         {2, 9},
         {1, 9},
         {1, 8},
         {1, 7}}};
    auto correct = correct_destinations.cbegin();
    for (const auto distance : distances) {
        generator.seed(42);
        int start_row = 8;
        int start_col = 7;
        if (!network.has_node_at(start_row, start_col)) {
            std::cerr << "Expected node at " << start_row << ", " << start_col << "\n";
            ret += 1;
            break;
        }
        int end_row;
        int end_col;
        std::tie(end_row, end_col) =
            network.travel(start_row, start_col, distance, generator);
        if (correct->first != end_row || correct->second != end_col) {
            std::cerr << "from (" << start_row << ", " << start_col << ") to ("
                      << end_row << ", " << end_col << ") in " << distance
                      << " but expected to arrive to (" << correct->first << ", "
                      << correct->second << ")\n";
            ret += 1;
        }
        ++correct;
    }
    if (ret) {
        std::cerr << "---\n";
        network.dump_yaml(std::cout);
    }
    return ret;
}

int test_snap_network()
{
    int ret = 0;
    BBox<double> bbox;
    bbox.north = 10;
    bbox.south = 0;
    bbox.east = 30;
    bbox.west = 20;
    Network<int> network{bbox, 1, 1, true};
    std::stringstream network_stream{
        "1,2,21.4;7.5;22.3;7.2\n"
        "1,4,21.4;7.5;21.9;8.0;22.5;8.6\n"
        "5,1,27.5;1.5;26.7;1.4;25.9;1.3;25.2;1.1;24.5;1.7;23.9;2.3;23.2;2.8;22.5;2.3;21.8;1.8;21.2;1.3;20.5;1.4;20.7;2.3;20.8;3.2;20.9;4.0;21.0;4.9;21.1;5.7;21.3;6.6;21.4;7.5\n"
        "2,8,22.3;7.2;23.2;7.1;24.0;6.9;24.8;6.7;25.7;6.6;26.5;6.4\n"
        "8,10,26.5;6.4;26.8;5.7;27.2;4.9;27.5;4.2;27.9;3.5;28.2;2.7\n"
        "11,8,28.3;9.0;28.5;8.1;28.6;7.3;28.2;6.8;27.7;6.3;27.1;6.4;26.5;6.4\n"
        "2,5,22.3;7.2;23.0;7.7;23.7;8.1;24.3;8.5;25.0;9.0;25.8;9.0;26.6;9.0;27.4;9.0;28.2;9.0;29.0;9.0;29.0;8.0;29.0;7.0;29.0;6.0;29.0;5.0;29.0;4.0;29.0;3.0;29.0;2.0;29.0;1.0;28.2;1.3;27.5;1.5\n"
        "5,8,27.5;1.5;26.7;1.8;26.0;2.0;26.1;2.9;26.2;3.8;26.3;4.7;26.4;5.5;26.5;6.4\n"};
    network.load(network_stream);

    std::default_random_engine generator;
    const int num_distances = 5;
    std::array<std::tuple<double, int, int>, num_distances> destinations{
        {{1, 1, 8}, {1.9, 1, 8}, {2.1, 3, 6}, {5.9, 3, 6}, {6.1, 2, 2}}};
    for (const auto destination : destinations) {
        generator.seed(42);
        int start_row = 1;
        int start_col = 8;
        if (!network.has_node_at(start_row, start_col)) {
            std::cerr << "Expected node at " << start_row << ", " << start_col << "\n";
            ret += 1;
            break;
        }
        int end_row;
        int end_col;
        std::tie(end_row, end_col) =
            network.travel(start_row, start_col, std::get<0>(destination), generator);
        if (std::get<1>(destination) != end_row
            || std::get<2>(destination) != end_col) {
            std::cerr << "from (" << start_row << ", " << start_col << ") to ("
                      << end_row << ", " << end_col << ") in "
                      << std::get<0>(destination) << " but expected to arrive to ("
                      << std::get<1>(destination) << ", " << std::get<2>(destination)
                      << ")\n";
            ret += 1;
        }
    }
    if (ret) {
        // network.dump_yaml(std::cout);
    }
    return ret;
}

int test_cost_network()
{
    int ret = 0;
    BBox<double> bbox;
    bbox.north = 10;
    bbox.south = 0;
    bbox.east = 30;
    bbox.west = 20;
    Network<int> network{bbox, 1, 1};
    std::stringstream network_stream{
        "node_1,node_2,cost,geometry\n"
        "1,2,111.5,21.4;7.5;22.3;7.2\n"
        "1,4,112.6,21.4;7.5;21.9;8.0;22.5;8.6\n"
        "2,8,214.6,22.3;7.2;23.2;7.1;24.0;6.9;24.8;6.7;25.7;6.6;26.5;6.4\n"
        "8,10,126.4,26.5;6.4;26.8;5.7;27.2;4.9;27.5;4.2;27.9;3.5;28.2;2.7\n"
        "11,8,127.4,28.3;9.0;28.5;8.1;28.6;7.3;28.2;6.8;27.7;6.3;27.1;6.4;26.5;6.4\n"
        "5,8,129.7,27.5;1.5;26.7;1.8;26.0;2.0;26.1;2.9;26.2;3.8;26.3;4.7;26.4;5.5;26.5;6.4\n"};
    network.load(network_stream);

    std::default_random_engine generator;
    const int num_distances = 9;
    std::array<std::tuple<double, int, int>, num_distances> destinations{
        {{70, 6, 6},
         {90, 5, 6},
         {110, 4, 6},
         {130, 3, 6},
         {150, 3, 6},
         {170, 3, 5},
         {190, 3, 5},
         {210, 3, 5},
         {230, 3, 4}}};
    for (const auto destination : destinations) {
        generator.seed(42);
        int start_row = 8;
        int start_col = 7;
        if (!network.has_node_at(start_row, start_col)) {
            std::cerr << "Expected node at " << start_row << ", " << start_col << "\n";
            ret += 1;
        }
        int end_row;
        int end_col;
        std::tie(end_row, end_col) =
            network.travel(start_row, start_col, std::get<0>(destination), generator);
        if (std::get<1>(destination) != end_row
            || std::get<2>(destination) != end_col) {
            std::cerr << "from (" << start_row << ", " << start_col << ") to ("
                      << end_row << ", " << end_col << ") in "
                      << std::get<0>(destination) << " but expected to arrive to ("
                      << std::get<1>(destination) << ", " << std::get<2>(destination)
                      << ")\n";
            ret += 1;
        }
    }
    if (ret) {
        // std::cout << "---\n";
        // network.dump_yaml(std::cout);
    }
    return ret;
}

int test_network_negative_probability()
{
    BBox<double> bbox;
    bbox.north = 10;
    bbox.south = 0;
    bbox.east = 30;
    bbox.west = 20;
    Network<int> network{bbox, 1, 1, true};
    std::stringstream network_stream{
        "node_1,node_2,probability,geometry"
        "1,2,0,21;7;22;7\n"
        "1,4,1,21;7;22;8\n"
        "5,1,-0.5,27;1;21;7\n"
        "2,8,-80,22;7;26;6\n"};
    try {
        network.load(network_stream);
        std::cerr << "Network loaded without an exception for negative probability\n";
        return 1;
    }
    catch (const std::invalid_argument& error) {
        // All is good. We expect and exception.
    }
    return 0;
}

int test_network_correct_probability()
{
    BBox<double> bbox;
    bbox.north = 10;
    bbox.south = 0;
    bbox.east = 30;
    bbox.west = 20;
    Network<int> network{bbox, 1, 1, true};
    std::stringstream network_stream{
        "node_1,node_2,probability,geometry"
        "1,2,0,21;7;22;7\n"
        "1,4,1,21;7;22;8\n"
        "5,1,0.4,27;1;21;7\n"
        "2,8,0.8,22;7;26;6\n"};
    try {
        network.load(network_stream);
        return 0;
    }
    catch (const std::invalid_argument& error) {
        std::cerr << "Network loaded with correct probability caused "
                     "an std::invalid_argument exception: "
                  << error.what() << "\n";
    }
    return 1;
}

int test_network_bad_probability_0_100()
{
    BBox<double> bbox;
    bbox.north = 10;
    bbox.south = 0;
    bbox.east = 30;
    bbox.west = 20;
    Network<int> network{bbox, 1, 1, true};
    std::stringstream network_stream{
        "node_1,node_2,probability,geometry"
        "1,2,0,21;7;22;7\n"
        "1,4,1,21;7;22;8\n"
        "5,1,5,27;1;21;7\n"
        "2,8,100,22;7;26;6\n"};
    try {
        network.load(network_stream);
        std::cerr << "Network loaded without an exception "
                     "for probability outside of [0,1]\n";
        return 1;
    }
    catch (const std::invalid_argument& error) {
        // All is good. We expect and exception.
    }
    return 0;
}

int test_step_network()
{
    int ret = 0;
    BBox<double> bbox;
    bbox.north = 10;
    bbox.south = 0;
    bbox.east = 30;
    bbox.west = 20;
    Network<int> network{bbox, 1, 1, true};
    std::stringstream network_stream{
        "1,2,21;7;22;7\n"
        "1,4,21;7;22;8\n"
        "5,1,27;1;21;7\n"
        "2,8,22;7;26;6\n"
        "8,9,26;6;28;2\n"
        "3,8,28;9;26;6\n"
        "2,5,22;7;27;1\n"
        "5,8,27;1;26;6\n"};
    network.load(network_stream);

    std::default_random_engine generator;
    const int num_distances = 5;
    std::array<std::tuple<int, int, int>, num_distances> destinations{
        {{1, 4, 6}, {2, 3, 2}, {3, 9, 7}, {4, 4, 6}, {5, 1, 8}}};
    for (const auto destination : destinations) {
        generator.seed(42);
        int start_row = 1;
        int start_col = 8;
        if (!network.has_node_at(start_row, start_col)) {
            std::cerr << "Expected node at " << start_row << ", " << start_col << "\n";
            ret += 1;
            break;
        }
        int end_row;
        int end_col;
        std::tie(end_row, end_col) =
            network.step(start_row, start_col, generator, std::get<0>(destination));
        if (std::get<1>(destination) != end_row
            || std::get<2>(destination) != end_col) {
            std::cerr << "from (" << start_row << ", " << start_col << ") to ("
                      << end_row << ", " << end_col << ") in "
                      << std::get<0>(destination) << " but expected to arrive to ("
                      << std::get<1>(destination) << ", " << std::get<2>(destination)
                      << ")\n";
            ret += 1;
        }
    }
    if (ret) {
        // network.dump_yaml(std::cout);
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
    Network<int> network{bbox, 0.01, 0.01};
    if (network.has_node_at(1, 1)) {
        std::cerr << "Empty network should not have a node at any cell\n";
        ret += 1;
    }
    std::stringstream network_stream{
        "1,2,-79.937;37.270;-79.936;37.270;-79.936;37.271;-79.936;37.271;-79.936;37.271;-79.934;37.272;-79.934;37.272\n"
        "3,4,-79.902;37.367;-79.903;37.366;-79.903;37.366;-79.904;37.366;-79.905;37.365;-79.905;37.36;-79.920;37.352;-79.93;37.273;-79.940;37.273;-79.941;37.273\n"};
    network.load(network_stream);
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

/** Placeholder for string to string no-op. */
std::string convert_to(const std::string& text, std::string tag)
{
    UNUSED(tag);
    return text;
}

/** Convert string to double */
double convert_to(const std::string& text, double tag)
{
    UNUSED(tag);
    return std::stod(text);
}

/** Convert string to int */
int convert_to(const std::string& text, int tag)
{
    UNUSED(tag);
    return std::stoi(text);
}

/** Convert string to bool */
int convert_to(const std::string& text, bool tag)
{
    UNUSED(tag);
    if (text == "true" || text == "True" || text == "TRUE")
        return true;
    if (text == "false" || text == "false" || text == "FALSE")
        return false;
    throw std::invalid_argument(
        std::string("Value conversion error: ") + text + " is not a boolean value");
}

/**
 * \brief A generic key-value configuration which can convert values to desired types.
 */
class RawConfig
{
public:
    /**
     * Get value from config or the default value.
     *
     * The type is determined by the template parameter is specified.
     */
    template<typename T = std::string>
    T get(const std::string& key) const
    {
        auto it{values_.find(key)};
        if (it != values_.end())
            return convert_to(it->second, T());
        throw std::invalid_argument(std::string("No value for key: ") + key);
    }
    /**
     * Get value from config or the default value.
     *
     * The type is determined by the default value unless the template parameter is
     * specified. For floating point numbers, you can do:
     *
     * ```
     * double x = config.get("x", 0.);
     * ```
     *
     * However, if `0` is passed, the function parse the value as an integer and return
     * an integer.
     */
    template<typename T = std::string>
    T get(const std::string& key, T default_value) const
    {
        auto it{values_.find(key)};
        if (it != values_.end())
            return convert_to(it->second, T());
        return default_value;
    }
    // std::optional for default_value can replace the get overload in C++17.
    template<typename T>
    void set(const std::string& key, const T& value)
    {
        values_[key] = value;
    }

protected:
    std::map<std::string, std::string> values_;
};

/**
 * \brief Create RawConfig from a YAML-like text stream of simple `key: value` pairs.
 */
template<typename Stream>
RawConfig read_config(Stream& stream)
{
    RawConfig config;
    std::string line;
    while (std::getline(stream, line)) {
        std::regex delimeter(R"([\s]*:[\s]*)");
        std::smatch match;
        if (regex_search(line, match, delimeter)) {
            config.set(match.prefix(), match.suffix());
        }
        else {
            throw std::runtime_error(std::string("Incorrect format at line: ") + line);
        }
    }
    return config;
}

/**
 * @brief Create bounding box from configuration
 *
 * @param configuration containing north, south, east, and west
 * @return bounding box object with doubles
 */
BBox<double> bbox_from_config(const RawConfig& config)
{
    BBox<double> bbox;
    try {
        bbox.north = config.get<double>("north");
        bbox.south = config.get<double>("south");
        bbox.east = config.get<double>("east");
        bbox.west = config.get<double>("west");
    }
    catch (std::out_of_range&) {
        throw std::runtime_error(
            "At least one of north, south, east, west keys is missing");
    }
    return bbox;
}

int create_network_from_files(int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0]
                  << " read|stats|write|trips|trace CONFIG_FILE NETWORK_FILE\n";
        return 1;
    }
    std::string command = argv[1];
    bool show_stats = false;
    bool write_network = false;
    bool trips = false;
    bool trace = false;
    if (command == "stats") {
        show_stats = true;
    }
    else if (command == "write") {
        write_network = true;
    }
    else if (command == "trips") {
        trips = true;
    }
    else if (command == "trace") {
        trace = true;
    }
    else if (command != "read") {
        std::cerr << "Unknown sub-command: " << command << "\n";
        std::cerr << "Supported sub-commands are: read, stats, write, trips, trace\n";
        return 1;
    }

    std::string config_file{argv[2]};
    std::ifstream config_stream{config_file};
    if (!config_stream.is_open()) {
        std::cerr << "Failed to open config file: " << config_file << "\n";
        return 1;
    }
    std::string network_file{argv[3]};
    std::ifstream network_stream{network_file};
    if (!network_stream.is_open()) {
        std::cerr << "Failed to open network file: " << network_file << "\n";
        return 1;
    }

    RawConfig config = read_config(config_stream);
    BBox<double> bbox = bbox_from_config(config);
    double nsres = config.get<double>("nsres");
    double ewres = config.get<double>("ewres");

    Network<int> network(bbox, nsres, ewres, config.get("snap", false));
    network.load(network_stream);

    if (show_stats) {
        auto stats = network.collect_statistics();
        for (const auto& item : stats) {
            // TODO: Resolve the stats output directly in the network.
            // Replace by starts_with in C++20.
            if (item.first.rfind("standalone_node_", 0) == 0) {
                int row;
                int col;
                std::tie(row, col) = network.get_node_row_col(item.second);
                std::cout << item.first << ":\n";
                std::cout << "  id: " << item.second << "\n";
                std::cout << "  row: " << row << "\n";
                std::cout << "  col: " << col << "\n";
            }
            else {
                std::cout << item.first << ": " << item.second << "\n";
            }
        }
    }
    if (write_network) {
        network.dump_yaml(std::cout);
    }
    if (trips || trace) {
        double min_distance = config.get("min_distance", 1.);
        double max_distance = config.get("max_distance", 1.);
        double distance_increment = config.get("distance_increment", 1.);
        int seed = config.get("seed", 1);
        std::default_random_engine generator;
        // Seed the generator for trips (only here) and for trace, seed it here
        // for node shuffle (or random selection when using a subset).
        generator.seed(seed);
        auto nodes = network.get_all_nodes();
        shuffle_container(nodes, generator);
        int num_nodes = config.get<int>("num_start_nodes", nodes.size());

        if (trace)
            std::cout << "traces:\n";
        else
            std::cout << "trips:\n";
        for (int i = 0; i < num_nodes; ++i) {
            const auto& node = nodes[i % nodes.size()];
            int start_row = node.second.first;
            int start_col = node.second.second;
            if (!network.has_node_at(start_row, start_col)) {
                std::cerr << "Internal error in the Network\n";
                return 1;
            }
            // For trace, the trips for node are the trace.
            // For trips and distance range, trips are independent combinations.
            std::vector<std::tuple<int, int, double>> trips;
            for (double distance = min_distance; distance <= max_distance;
                 distance += distance_increment) {
                if (trace)
                    generator.seed(seed);
                int end_row;
                int end_col;
                std::tie(end_row, end_col) =
                    network.travel(start_row, start_col, distance, generator);
                trips.emplace_back(end_row, end_col, distance);
            }
            if (trace) {
                std::cout << "  - cells: [";
                for (const auto& cell : trips) {
                    std::cout << "[";
                    std::cout << std::get<0>(cell) << ", " << std::get<1>(cell);
                    std::cout << "], ";
                }
                std::cout << "]\n";
                std::cout << "    distances: [";
                for (const auto& cell : trips) {
                    std::cout << std::get<2>(cell) << ", ";
                }
                std::cout << "]\n";
                std::cout << "    start:\n";
                std::cout << "      row: " << start_row << "\n";
                std::cout << "      col: " << start_col << "\n";
                std::cout << "    end:\n";
                std::cout << "      row: " << std::get<0>(trips.back()) << "\n";
                std::cout << "      col: " << std::get<1>(trips.back()) << "\n";
            }
            else {
                for (const auto& cell : trips) {
                    std::cout << "  - start:\n";
                    std::cout << "      row: " << start_row << "\n";
                    std::cout << "      col: " << start_col << "\n";
                    std::cout << "    end:\n";
                    std::cout << "      row: " << std::get<0>(cell) << "\n";
                    std::cout << "      col: " << std::get<1>(cell) << "\n";
                    std::cout << "    distance: " << std::get<2>(cell) << "\n";
                }
            }
        }
    }

    return 0;
}

int run_tests()
{
    int ret = 0;

    ret += test_bbox_functions();
    ret += test_create_network();
    ret += test_travel_network();
    ret += test_snap_network();
    ret += test_cost_network();
    ret += test_step_network();
    ret += test_network_bad_probability_0_100();
    ret += test_network_negative_probability();
    ret += test_network_correct_probability();

    if (ret)
        std::cerr << "Number of errors in the network test: " << ret << "\n";
    return ret;
}

int main(int argc, char** argv)
{
    if (argc > 1)
        return create_network_from_files(argc, argv);
    else
        return run_tests();
}
