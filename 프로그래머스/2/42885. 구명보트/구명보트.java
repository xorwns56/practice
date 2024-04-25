import java.util.Arrays;
class Solution {
    public int solution(int[] people, int limit){
        Arrays.sort(people);
        int n = people.length - 1;
        int count = 0;
        for(int i=0;i<=n;){
            if(people[i]+people[n]>limit) n--;
            else{
                i++;
                n--;
            }
            count++;
        }
        return count;
    }
}