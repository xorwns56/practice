class Solution {
    public int solution(int[][] sizes) {
        int total_v = 0;
        int total_h = 0;
        for(int i = 0; i < sizes.length; i++){
            int v = Math.min(sizes[i][0], sizes[i][1]);
            int h = Math.max(sizes[i][0], sizes[i][1]);
            total_v = Math.max(total_v, v);
            total_h = Math.max(total_h, h);
        }
        return total_v * total_h;
    }
}