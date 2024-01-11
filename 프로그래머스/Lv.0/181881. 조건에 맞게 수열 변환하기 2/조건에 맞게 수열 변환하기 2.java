class Solution {
    public int solution(int[] arr) {
        int iter = 0;
        while(true){
            int count = 0;
            for(int i = 0; i < arr.length; i++){
                if(arr[i] >= 50 && (arr[i] & 1) == 0) arr[i] /= 2;
                else if(arr[i] < 50 && (arr[i] & 1) == 1) arr[i] = arr[i] * 2 + 1;
                else count++;
            }
            if(count == arr.length) return iter;
            iter++;
        }
    }
}