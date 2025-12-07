#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

//hello hanyong

//Parking space dimensions
//consider bottom left corner to be (0,0)
const double SPACE_WIDTH  = 2.5;   // horizontal (x direction)
const double SPACE_LENGTH = 5.0;  //vertical (y direction)

//define constants
const double ONE_WAY_ENTRANCE = 4.5;
const double TWO_WAY_ENTRANCE = 7.0;
const double BUFFER = 2.5;
const double ONE_WAY_LANE = 4.5;

//minimum lot size
const double MIN_WIDTH  = 14.5;
const double MIN_LENGTH = 7.5;


//Get user dimensions
void getDimensions(double &width, double &length) {
    cout << "Enter parking lot width (x direction): ";
    cin >> width;

    cout << "Enter parking lot length (y direction): ";
    cin >> length;

    cout << fixed << setprecision(2);
}



//Check entrance size
int checkEntr(double width) {
    if (width < ONE_WAY_ENTRANCE) {
        return 0;  //no entrance possible
    }
    if (width >= TWO_WAY_ENTRANCE) {
        return 2;  //two way available
    }
    return 1;      //only one way
}



//Check minimum size
bool checkMinSize(double width, double length) {
    return (width >= MIN_WIDTH && length >= MIN_LENGTH);
}


//Generate parking spaces
//Return vector of (x, y) pairs
vector<pair<double,double>> generateParkingSpaces(double lotWidth) {
    vector<pair<double,double>> coords;

    int spacesPerLane = (int)(lotWidth / SPACE_WIDTH);
    int numberOfLanes = 2;

    double firstLaneY = BUFFER;
    double secondLaneY = firstLaneY + SPACE_LENGTH + ONE_WAY_LANE;

    for (int lane = 0; lane < numberOfLanes; lane++) {
        double y_base = firstLaneY + lane * (SPACE_LENGTH + ONE_WAY_LANE);

        for (int i = 0; i < spacesPerLane; i++) {
            double x_right = (i + 1) * SPACE_WIDTH;
            double y_right = y_base;

            coords.push_back({x_right, y_right});
        }
    }

    return coords;
}


//main program
int main() {
    double width, length;

    //input
    getDimensions(width, length);

    //check entrance type
    int entranceType = checkEntr(width);

    if (entranceType == 0) {
        cout << "The width too small for a one-way entrance.\n";
        return 0;
    }
    else if (entranceType == 2) {
        cout << "Two way entrance possible (" << TWO_WAY_ENTRANCE << " m).\n";
    }
    else {
        cout << "Only one-way entrance possible (" << ONE_WAY_ENTRANCE << " m).\n";
    }

    //Minimum lot size check
    if (!checkMinSize(width, length)) {
        cout << "Lot size is too small. A parking lot is defined to be at least a one way entrance and 2 lanes of parking spaces seperated by a one way lane.\n";
        return 0;
    }

    //Generate spaces
    vector<pair<double,double>> spaces = generateParkingSpaces(width);

    //Output results
    cout << "Bottom-right coordinates of each parking space:\n";
    for (size_t i = 0; i < spaces.size(); i++) {
        cout << "(" << spaces[i].first << ", " << spaces[i].second << ")\n";
    }

    return 0;
}
