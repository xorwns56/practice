import java.util.*;
class Bridge{
    int time;
    int weight;
    int length;
    Queue<Truck> trucks;
    Queue<Integer> truck_weights;
    int curr_weight;
    Bridge(int length, int weight, int[] truck_weights){
        this.time = 0;
        this.length = length;
        this.weight = weight;
        this.curr_weight = 0;
        this.trucks = new LinkedList<>();
        this.truck_weights = new LinkedList<>();
        for(int i = 0; i < truck_weights.length; i++) this.truck_weights.add(truck_weights[i]);
    }
    void tick(){
        time++;
        if(!trucks.isEmpty() && trucks.peek().out_time == time) curr_weight -= trucks.remove().weight;
        if(!truck_weights.isEmpty() && curr_weight + truck_weights.peek() <= weight){
            int truck_weight = truck_weights.remove();
            trucks.add(new Truck(truck_weight, time + length));
            curr_weight += truck_weight;
        }
    }
}

class Truck{
    int weight;
    int out_time;
    Truck(int weight, int out_time){
        this.weight = weight;
        this.out_time = out_time;
    }
}

class Solution {
    public int solution(int bridge_length, int weight, int[] truck_weights) {
        Bridge bridge = new Bridge(bridge_length, weight, truck_weights);
        while(!bridge.truck_weights.isEmpty() || !bridge.trucks.isEmpty()) bridge.tick();
        return bridge.time;
    }
}