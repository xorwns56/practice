import java.util.*;
class Solution {
    public int solution(int[][] dots) {
        Arrays.sort(dots, (x, y)-> x[0] == y[0] ? x[1] - y[1] : x[0] - y[0]);
        return (dots[2][0] - dots[0][0]) * (dots[1][1] - dots[0][1]);
    }
}