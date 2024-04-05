class Solution {
    public int solution(int[][] triangle) {
        int[][] sum = new int[triangle.length][triangle.length];
        return visit(triangle, sum, 0, 0);
    }
    public int visit(int[][] triangle, int[][] sum, int level, int index){
        if(level == triangle.length) return 0;
        if(sum[level][index] > 0) return sum[level][index];
        return sum[level][index] = triangle[level][index] + Math.max(visit(triangle, sum, level + 1, index), visit(triangle, sum, level + 1, index + 1));
    }
}