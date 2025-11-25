#include <iostream>
#include <vector>
#include <iomanip>
#include <utility>
#include <cmath>

using namespace std;

//dear hanyong omar jihee aethan
//kys

//why did the chicken cross the road
//cos it wanted to get to the other side
//haha

//(0,0) is considered to be bottom left corner of parking lot
//Constants 
const double SPACE_W = 2.5;   //width is x direction
const double SPACE_L = 5.0;   //length is y-direction

const double ONE_WAY_TRAFFIC = 4.5; //one way traffic lane dim
const double TWO_WAY_TRAFFIC = 7.0; //two way traffic lane dim

//min lot dimensions for initial test
//parking lot is defined as at least 1 way entrance/exit and two lanes of at least two parking spaces with a one way traffic lane
const double MIN_WIDTH  = 14.5;
const double MIN_LENGTH = 7.5;


//i got this from chat 
//realised it would be faster this way
struct Space {
    // bottom-right and top-left corners
    pair<double,double> bottomRight; // (x, y)
    pair<double,double> topLeft;     // (x, y)
};


//void function to get dimensions from user
void getDimensions(double &width, double &length) {
    cout << "Enter parking lot width (x direction, metres): ";
    cin  >> width;
    cout << "Enter parking lot length (y direction, metres): ";
    cin  >> length;
    cout << fixed << setprecision(3);
}

//compute maximum columns of parking spaces in lot
int maxParkingColumns(double lotWidth, double trafficW = ONE_WAY_TRAFFIC) {
    //maximum possible columns if no vertical traffic lanes
    int maxPossible = (int)floor(lotWidth / SPACE_W);
    //but there are so yadadadadada
    for (int C = maxPossible; C >= 1; --C) {
        double used = C * SPACE_W + (C - 1) * trafficW;
        if (used <= lotWidth + 1e-9) return C;
    }
    return 0;
}



//find max number of rows that can fit
//must have top and bottom lanes for connecting traffic
int maxParkingRows(double lotLength, double horizontalTraffic = ONE_WAY_TRAFFIC) {
    //available length after top and bottom becomes traffic lanes
    if (lotLength < 2*ONE_WAY_TRAFFIC + SPACE_L - 1e-9) return 0; // not enough
    //max rows
    int maxPossible = (int)floor( (lotLength - 2*ONE_WAY_TRAFFIC) / SPACE_L );
    for (int R = maxPossible; R >= 1; --R) {
        double used = 2*ONE_WAY_TRAFFIC + R*SPACE_L + (R - 1) * horizontalTraffic;
        if (used <= lotLength + 1e-9) return R;
    }
    return 0;
}


//generates coordinates for bottom right and top left coordinates
//chat gave me directions cus i am too lazy to do without rn
vector<Space> generateSpaces(double lotWidth, double lotLength,
                             int &outCols, int &outRows,
                             double trafficVertical = ONE_WAY_TRAFFIC,
                             double trafficHorizontal = ONE_WAY_TRAFFIC) {
    vector<Space> result;

    //find maximum columns and rows
    int C = maxParkingColumns(lotWidth, trafficVertical);
    int R = maxParkingRows(lotLength, trafficHorizontal);

    outCols = C;
    outRows = R;
    
    //just to check if theres no room for columns or rows
    if (C <= 0 || R <= 0) return result;

    //this fuckery is for finding the x and y coordinates
    //chat obviously
    //first column is along x = 0
    for (int i = 0; i < C; ++i) {
        double x_left = i * (SPACE_W + trafficVertical);
        double x_right = x_left + SPACE_W;

        //find y coordinates
        for (int j = 0; j < R; ++j) {
            double y_bottom = ONE_WAY_TRAFFIC + j * (SPACE_L + trafficHorizontal);
            double y_top = y_bottom + SPACE_L;

            Space s;
            s.bottomRight = { x_right, y_bottom };
            s.topLeft     = { x_left,  y_top   };
            result.push_back(s);
        }
    }

    return result;
}

//main function
int main() {
    double lotW, lotL;
    getDimensions(lotW, lotL);

    //checks whether width is too small, as the thingy will fail
    //this if checks whether theres room for at least a one way traffic lane
    if (lotW < ONE_WAY_TRAFFIC - 1e-9) {
        cout << "Error: width too small for traffic lanes\n";
        return 0;
    }
    //ts checks whether the inputted length and width are bigger than the minimum constant
    if (lotW < MIN_WIDTH || lotL < MIN_LENGTH) {
        cout << "Error: input smaller than the minimum dimensions\n";
        
    }

    //ts was chatted cus i am too lazy to do it my slef
    //it is currently 2 am i will kms
    //omar ik u reading this
    // Compute optimal grid packing using the minimal (one-way) traffic lane widths
    int cols = 0, rows = 0;
    //basically just finds columns 
    vector<Space> spaces = generateSpaces(lotW, lotL, cols, rows,
                                          ONE_WAY_TRAFFIC, ONE_WAY_TRAFFIC);


    //outputs total columns, rows, and spaces
    cout << "  Parking columns (across width): " << cols << "\n";
    cout << "  Parking rows    (along length): " << rows << "\n";
    cout << "  Total spaces: " << (int)spaces.size();

    //outputs bottom right and top left corner coordinates of each space
    cout << "List of parking spaces (bottom right and top left):\n";
    for (size_t i = 0; i < spaces.size(); ++i) {
        cout << "BR(" << spaces[i].bottomRight.first  << ", " << spaces[i].bottomRight.second << "), "
             << "TL(" << spaces[i].topLeft.first << ", "  << spaces[i].topLeft.second << ")\n";
    }

    //output used and unused space
    double usedWidth = cols * SPACE_W + (cols - 1) * ONE_WAY_TRAFFIC;
    double usedLength = 2*ONE_WAY_TRAFFIC + rows * SPACE_L + (rows - 1) * ONE_WAY_TRAFFIC;
    cout << "Space used (approx): width = " << usedWidth << " m, length = " << usedLength << " m\n";
    cout << "Unused space: width unused = " << max(0.0, lotW - usedWidth)
         << " m, length unused = " << max(0.0, lotL - usedLength) << " m\n";

    return 0;
}
