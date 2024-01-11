class Solution {
    public int solution(int[] num_list) {
        int count = 0;
        for(int i = 0; i < num_list.length; i++){
            while(num_list[i] > 1){
                count++;
                num_list[i] >>= 1;
            }
        }
        return count;
    }
}