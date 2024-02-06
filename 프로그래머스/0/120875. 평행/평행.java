
class Solution {
    public int solution(int[][] dots) {
        double x = (double)Math.abs(dots[0][0] - dots[1][0]) / Math.abs(dots[0][1] - dots[1][1]);
        double y = (double)Math.abs(dots[2][0] - dots[3][0]) / Math.abs(dots[2][1] - dots[3][1]);
        if(x == y) return 1;
        x = (double)Math.abs(dots[0][0] - dots[2][0]) / Math.abs(dots[0][1] - dots[2][1]);
        y = (double)Math.abs(dots[1][0] - dots[3][0]) / Math.abs(dots[1][1] - dots[3][1]);
        if(x == y) return 1;
        x = (double)Math.abs(dots[0][0] - dots[3][0]) / Math.abs(dots[0][1] - dots[3][1]);
        y = (double)Math.abs(dots[1][0] - dots[2][0]) / Math.abs(dots[1][1] - dots[2][1]);
        if(x == y) return 1;
        return 0;
    }
    
    
}