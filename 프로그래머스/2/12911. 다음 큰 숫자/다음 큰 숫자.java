class Solution {
    public int solution(int n) {
        int count = 0;
        int tmp = n;
        while(tmp > 0){
            if(tmp % 2 == 1) count++;
            tmp /= 2;
        }
        int x = 0;
        while(true){
            tmp = n + ++x;
            int tmp_count = 0;
            while(tmp > 0){
                if(tmp % 2 == 1) tmp_count++;
                tmp /= 2;
            }
            if(tmp_count == count) return n + x;
        }
    }
}