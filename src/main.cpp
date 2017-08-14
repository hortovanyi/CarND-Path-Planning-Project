#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() {
  return M_PI;
}
double deg2rad(double x) {
  return x * pi() / 180;
}
double rad2deg(double x) {
  return x * 180 / pi();
}

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x,
                    vector<double> maps_y) {

  double closestLen = 100000;  //large number
  int closestWaypoint = 0;

  for (int i = 0; i < maps_x.size(); i++) {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x, y, map_x, map_y);
    if (dist < closestLen) {
      closestLen = dist;
      closestWaypoint = i;
    }

  }

  return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x,
                 vector<double> maps_y) {

  int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2((map_y - y), (map_x - x));

  double angle = abs(theta - heading);

  if (angle > pi() / 4) {
    closestWaypoint++;
  }

  return closestWaypoint;

}



// get the lane based on frenet coordinates
int LaneFrenet(double d) {
  int lane_width=4;
  return int(floor(d/lane_width));
}

double FrenetLaneCenter(int lane) {
  int lane_width=4;
  return double((lane)*lane_width + lane_width/2);
}

double Normalise(double x) {
  return 2.0f / (1.0f + exp(-x)) - 1.0f;
}

vector<int> SensorFussionLaneIds(int lane, json sensor_fusion) {
  vector<int> ids;
  for (int i = 0; i < sensor_fusion.size(); i++) {
    float other_d = sensor_fusion[i][6];

    // see if the vehicle is in my lane
    int other_lane = LaneFrenet(other_d);
    if (other_lane < 0 || other_lane > 2)
      continue;

    if (other_lane == lane)
      ids.push_back(i);
  }
  return ids;
}

double NearestApproach(vector<int> ids, json sensor_fusion, double check_dist, double car_s) {
  double closest = 99999;

  for (int id : ids) {
    double vx = sensor_fusion[id][3];
    double vy = sensor_fusion[id][4];
    double check_speed = sqrt(vx * vx + vy * vy);
    double check_start_car_s = sensor_fusion[id][5];

    double check_end_car_s = check_start_car_s + check_dist * check_speed;

    double dist_start = fabs(check_start_car_s-car_s);
    if ( dist_start < closest)
      closest = dist_start;

    double dist_end = fabs(check_end_car_s-car_s);
    if ( dist_end < closest)
      closest = dist_end;
  }
  return closest;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta,
                         vector<double> maps_x, vector<double> maps_y) {
  int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

  int prev_wp;
  prev_wp = next_wp - 1;
  if (next_wp == 0) {
    prev_wp = maps_x.size() - 1;
  }

  double n_x = maps_x[next_wp] - maps_x[prev_wp];
  double n_y = maps_y[next_wp] - maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
  double proj_x = proj_norm * n_x;
  double proj_y = proj_norm * n_y;

  double frenet_d = distance(x_x, x_y, proj_x, proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000 - maps_x[prev_wp];
  double center_y = 2000 - maps_y[prev_wp];
  double centerToPos = distance(center_x, center_y, x_x, x_y);
  double centerToRef = distance(center_x, center_y, proj_x, proj_y);

  if (centerToPos <= centerToRef) {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for (int i = 0; i < prev_wp; i++) {
    frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
  }

  frenet_s += distance(0, 0, proj_x, proj_y);

  return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s,
                     vector<double> maps_x, vector<double> maps_y) {
  int prev_wp = -1;

  while (s > maps_s[prev_wp + 1] && (prev_wp < (int) (maps_s.size() - 1))) {
    prev_wp++;
  }

  int wp2 = (prev_wp + 1) % maps_x.size();

  double heading = atan2((maps_y[wp2] - maps_y[prev_wp]),
                         (maps_x[wp2] - maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s - maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
  double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

  double perp_heading = heading - pi() / 2;

  double x = seg_x + d * cos(perp_heading);
  double y = seg_y + d * sin(perp_heading);

  return {x,y};

}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from

  string map_file_ = "../data/highway_map.csv";
//  string map_file_ = "../../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }

  cout << "loaded " << map_waypoints_x.size() << " waypoints" << endl;

  // start in lane 1
  int lane = 1;

  // have a reference velocity to target
  double ref_vel = 0; // mph
  double max_vel = 49.75;

  h.onMessage(
      [&lane, &ref_vel, &max_vel, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy]
       (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
          uWS::OpCode opCode) {
        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        //auto sdata = string(data).substr(0, length);
        //cout << sdata << endl;
        if (length && length > 2 && data[0] == '4' && data[1] == '2') {

          auto s = hasData(data);

          if (s != "") {
            auto j = json::parse(s);

            string event = j[0].get<string>();

            if (event == "telemetry") {
              // j[1] is the data JSON object

              // Main car's localization Data
              double car_x = j[1]["x"];
              double car_y = j[1]["y"];
              double car_s = j[1]["s"];
              double car_d = j[1]["d"];
              double car_yaw = j[1]["yaw"];
              double car_speed = j[1]["speed"];

              cout << "json x " << car_x << " y "<< car_y << " s " << car_s << " d " << car_d << " yaw " << car_yaw << " speed " << car_speed << endl;

              // Previous path data given to the Planner
              auto previous_path_x = j[1]["previous_path_x"];
              auto previous_path_y = j[1]["previous_path_y"];
              // Previous path's end s and d values
              double end_path_s = j[1]["end_path_s"];
              double end_path_d = j[1]["end_path_d"];

              // Sensor Fusion Data, a list of all other cars on the same side of the road.
              auto sensor_fusion = j[1]["sensor_fusion"];

              int prev_size = previous_path_x.size();

              if (prev_size>0) {
                car_s = end_path_s;
              }


              // calculate lane speeds
              vector<double> lane_speeds ={0.0,0.0,0.0};
              vector<int> lane_count={0,0,0};

              for (int i = 0; i < sensor_fusion.size(); i++) {
                float other_d = sensor_fusion[i][6];
                int other_lane = LaneFrenet(other_d);
//                cout << " other d " << other_d << "other_lane " << other_lane;
                if (other_lane < 0 || other_lane > 2)
                   continue;
                assert(other_lane >= 0 && other_lane <=2 );
                double vx = sensor_fusion[i][3];
                double vy = sensor_fusion[i][4];
                double check_speed = sqrt(vx*vx+vy*vy);
                lane_speeds[other_lane]+=check_speed*2.23694; //convert to MPH
                lane_count[other_lane]+=1;
              }

              // lets average now
              for (int i = 0; i < lane_speeds.size(); i++) {
                int n = lane_count[i];
                // if no other vehicles in the lane - set to maximum allowed
                if (n == 0){
                  lane_speeds[i]=max_vel;
                }
                // average the speed and
                else {
                  lane_speeds[i]=lane_speeds[i]/n;
                }
              }

              bool too_close = false;
              double closest_speed = max_vel;

              for (int i = 0; i < sensor_fusion.size(); i++) {
                float other_d = sensor_fusion[i][6];

                // see if the vehicle is in my lane
                int other_lane = LaneFrenet(other_d);
                if (other_lane < 0 || other_lane > 2)
                  continue;

                int ego_lane = LaneFrenet(end_path_d);
//                cout << "ego d " << end_path_d << " lane "<< ego_lane << " other d " << other_d << " lane " << other_lane << endl;
                if (lane == other_lane ) {
                  double vx = sensor_fusion[i][3];
                  double vy = sensor_fusion[i][4];
                  double check_speed = sqrt(vx*vx+vy*vy);
                  double check_car_s = sensor_fusion[i][5];

                  check_car_s+=((double)prev_size*.02*check_speed);
                  // check s values greater (its ahead) and less than 30 gap
                  if (check_car_s > car_s && check_car_s-car_s < 30)
                  {
                    // TODO work out what speed or turn left or turn right
                    too_close = true;
                    closest_speed = check_speed;
                  }
                }
              }

              if (too_close) {
                vector<int> possible_lanes;
                if (lane == 0)
                  possible_lanes = {0,1};
                else if (lane == 1)
                  possible_lanes ={0,1,2};
                else
                  possible_lanes ={1,2};

                int best_lane=lane;
                double best_cost=numeric_limits<double>::max();
                for (int check_lane: possible_lanes) {
                  double cost = 0;

                  // lane cost
                  if (check_lane != lane)
                    cost += 1000;

                  // speed cost
                  double avg_speed =lane_speeds[check_lane];
                  cost += Normalise(2.0 * (avg_speed-ref_vel/avg_speed)) * 1000;

                  // collision and buffer cost
                  auto ids = SensorFussionLaneIds(check_lane,sensor_fusion);
                  double nearest = NearestApproach(ids, sensor_fusion, .02*prev_size, car_s);

                  double buffer = 10;

                  if (nearest < buffer)
                    cost+=pow(10,5);

                  cost += Normalise(2*buffer/nearest) * 1000;

                  if (cost < best_cost) {
                    best_lane = check_lane;
                    best_cost = cost;
                  }
                }

                // if we are close and going faster then the average
                if (best_lane == lane && (ref_vel > lane_speeds[lane] || ref_vel > closest_speed))
                  ref_vel -= .224;

                // change lanes
                lane = best_lane;

              } else
              if (ref_vel < max_vel)  // MPH
              {
                ref_vel += .224;
              }

              vector<double> ptsx, ptsy;

              // reference x,y, yaw states
              double ref_x = car_x;
              double ref_y = car_y;
              double ref_yaw = deg2rad(car_yaw);

              // if previous size is almost empty, use the car as starting reference
              if(prev_size<2)
              {
                // Use the two points the make the path tangent to the car
                double prev_car_x = car_x - cos(car_yaw);
                double prev_car_y = car_y - sin(car_yaw);

                ptsx.push_back(prev_car_x);
                ptsx.push_back(car_x);

                ptsy.push_back(prev_car_y);
                ptsy.push_back(car_y);
              }
              // use the previous paths end point as a starting reference
              else
              {
                ref_x = previous_path_x[prev_size-1];
                ref_y = previous_path_y[prev_size-1];

                double ref_x_prev = previous_path_x[prev_size-2];
                double ref_y_prev = previous_path_y[prev_size-2];

                ref_yaw = atan2(ref_y-ref_y_prev,ref_x-ref_x_prev);

                ptsx.push_back(ref_x_prev);
                ptsx.push_back(ref_x);

                ptsy.push_back(ref_y_prev);
                ptsy.push_back(ref_y);
              }

//              int lane_d = 2+4*lane;
              int lane_d = FrenetLaneCenter(lane);
              vector<double> next_wp0 = getXY(car_s+30,lane_d,map_waypoints_s, map_waypoints_x, map_waypoints_y );
              vector<double> next_wp1 = getXY(car_s+60,lane_d,map_waypoints_s, map_waypoints_x, map_waypoints_y );
              vector<double> next_wp2 = getXY(car_s+90,lane_d,map_waypoints_s, map_waypoints_x, map_waypoints_y );

              ptsx.push_back(next_wp0[0]);
              ptsx.push_back(next_wp1[0]);
              ptsx.push_back(next_wp2[0]);

              ptsy.push_back(next_wp0[1]);
              ptsy.push_back(next_wp1[1]);
              ptsy.push_back(next_wp2[1]);

//              cout << "global pts x ";
//              for (auto n: ptsx) cout << n << " ";
//              cout << endl;
//
//              cout << "global pts y ";
//              for (auto n: ptsy) cout << n << " ";
//              cout << endl;

              for (int i = 0; i < ptsx.size(); i++)
              {
                // shift car reference angle to 0 degrees
                double shift_x = ptsx[i]-ref_x;
                double shift_y = ptsy[i]-ref_y;

                ptsx[i] = (shift_x * cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));
                ptsy[i] = (shift_x * sin(0-ref_yaw)+shift_y*cos(0-ref_yaw));
              }

//              cout << "car pts x ";
//              for (auto n: ptsx) cout << n << " ";
//              cout << endl;
//
//              cout << "car pts y ";
//              for (auto n: ptsy) cout << n << " ";
//              cout << endl;

              // create a spline
              tk::spline s;

              // set (x,y) points to the spline
              s.set_points(ptsx,ptsy);

              vector<double> next_x_vals, next_y_vals;
              // start with all the previous path points from the last time
              for(int i = 0; i < previous_path_x.size(); i++) {
                next_x_vals.push_back(previous_path_x[i]);
                next_y_vals.push_back(previous_path_y[i]);
              }

              // break the spline up
              double target_x = 30.0;
              double target_y = s(target_x);
              double target_dist = sqrt(pow(target_x,2)+pow(target_y,2));

              double x_add_on = 0;

              for(int i = 1; i < 50-previous_path_x.size(); i++) {
                double N = (target_dist/(.02*ref_vel/2.24));
                double x_point = x_add_on+(target_x)/N;
                double y_point = s(x_point);

                x_add_on = x_point;

                double x_ref = x_point;
                double y_ref = y_point;

                // rotate back to normal after rotating to vehicle space
                x_point = (x_ref * cos(ref_yaw)-y_ref*sin(ref_yaw));
                y_point = (x_ref * sin(ref_yaw)+y_ref*cos(ref_yaw));

                x_point += ref_x;
                y_point += ref_y;

                next_x_vals.push_back(x_point);
                next_y_vals.push_back(y_point);
              }

//              cout << "next x ";
//              for (auto n: next_x_vals) cout << n << " ";
//              cout << endl;
//
//              cout << "next y ";
//              for (auto n: next_y_vals) cout << n << " ";
//              cout << endl;
//              double dist_inc = 0.3;
//              for(int i = 0; i < 50; i++)
//              {
//                double next_s = car_s+(i+1)*dist_inc;
//                double next_d = 6;
//                vector<double> xy = getXY(next_s, next_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
////                cout << "next_s " << next_s << " next_d " << next_d << " x " << xy[0] << " y " << xy[1] << endl;
//                next_x_vals.push_back(xy[0]);
//                next_y_vals.push_back(xy[1]);
//              }

              json msgJson;
              msgJson["next_x"] = next_x_vals;
              msgJson["next_y"] = next_y_vals;

              auto msg = "42[\"control\","+ msgJson.dump()+"]";

              //this_thread::sleep_for(chrono::milliseconds(1000));
              ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

            }
          } else {
            // Manual driving
            std::string msg = "42[\"manual\",{}]";
            ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          }
        }
      });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
      size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
      char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}

