class Solution {
    public int solution(int[] array, int n) {
        int min_dist_idx = -1;
        int min_dist = Integer.MAX_VALUE;
        for(int i = 0; i < array.length; i++){
            int dist = Math.abs(array[i] - n);
            if(min_dist > dist){
                min_dist_idx = i;
                min_dist = dist;
            }else if(min_dist == dist && array[min_dist_idx] > array[i]) min_dist_idx = i;
        }
        return array[min_dist_idx];
    }
}