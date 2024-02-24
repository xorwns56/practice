import java.util.*;
class Solution {
    public int[] solution(String[] park, String[] routes) {
        int[] pos = new int[] { 0, 0 };
        boolean[][] arr = new boolean[park.length][park[0].length()];
        for(int i = 0; i < park.length; i++){
            for(int j = 0; j < park[i].length(); j++){
                char c = park[i].charAt(j);
                if(c == 'S'){
                    pos[0] = i;
                    pos[1] = j;
                }else if(c == 'X') arr[i][j] = true;
            }
        }
        for(int i = 0; i < routes.length; i++){
            char op = routes[i].charAt(0);
            int n = routes[i].charAt(2) - '0';
            int x = pos[1];
            int y = pos[0];
            boolean move_fail = false;
            for(int move = 0; move < n; move++){
                switch(op){
                    case 'N': y--; break;
                    case 'E': x++; break;
                    case 'S': y++; break;
                    case 'W': x--; break;
                }
                if(!(0 <= y && y < arr.length && 0 <= x && x < arr[0].length) || arr[y][x]){
                    move_fail = true;
                    break;
                }
            }
            if(!move_fail){
                pos[0] = y;
                pos[1] = x;
            }
        }
        return pos;
    }
}