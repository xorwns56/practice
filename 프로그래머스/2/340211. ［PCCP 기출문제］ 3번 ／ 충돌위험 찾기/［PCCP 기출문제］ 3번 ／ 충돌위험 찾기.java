import java.util.*;
class Solution {
    public int solution(int[][] points, int[][] routes) {
        int answer = 0;
        List<List<int[]>> routeList = new ArrayList<>();
        for(int i = 0; i < routes.length; i++){
            List<int[]> route = new ArrayList<>();
            for(int j = 0; j < routes[i].length; j++){
                route.add(new int[]{points[routes[i][j] - 1][0], points[routes[i][j] - 1][1]});
            }
            routeList.add(route);
        }
        int[][] curr = new int[routeList.size()][2];
        int[][] currCount = new int[200][200];
        for(int i = 0; i < routeList.size(); i++){
            curr[i] = routeList.get(i).remove(0);
            currCount[curr[i][0]][curr[i][1]]++;
            if(currCount[curr[i][0]][curr[i][1]] == 2){
                answer++;
            }
        }
        int completeCount = 0;
        while(completeCount < routeList.size()){
            currCount = new int[200][200];
            for(int i = 0; i < routeList.size(); i++){
                List<int[]> route = routeList.get(i);
                if(route.size() > 0){
                    int[] currPos = curr[i];
                    int[] destPos = route.get(0);
                    if(currPos[0] != destPos[0]){
                        currPos[0] += currPos[0] < destPos[0] ? 1 : -1; 
                    }else{
                        currPos[1] += currPos[1] < destPos[1] ? 1 : -1;
                    }
                    currCount[currPos[0]][currPos[1]]++;
                    if(currCount[currPos[0]][currPos[1]] == 2){
                        answer++;
                    }
                    if(currPos[0] == destPos[0] && currPos[1] == destPos[1]){
                        route.remove(0);
                        if(route.size() == 0) completeCount++;
                    }
                }
            }
        }
        return answer;
    }
}