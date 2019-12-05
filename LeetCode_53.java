/*
53. Maxinum Subarray
Given an integer array nums, find the contiguous subarray (containing at least one number) which has the largest sum and return its sum.

Example:
Input: [-2,1,-3,4,-1,2,1,-5,4],
Output: 6
Explanation: [4,-1,2,1] has the largest sum = 6.
*/

int maxSubArray(int* nums, int numsSize){
    int table[numsSize];
    int subArraySize = 1; //SubArray의 크기는 최소 1이므로 1로 초기화
    int max = -2147483647; //max를 int형의 최솟값으로 초기화
    
    while(subArraySize<=numsSize){ //SubArray의 크기가 numsSize 이하일 동안 반복
        for(int i=0;i<=numsSize-subArraySize;i++){
            if(subArraySize==1) table[i] = nums[i];
            else table[i] += nums[i+subArraySize-1];
            if(max < table[i]) max = table[i]; //max값 갱신
        }
        subArraySize++;
    }
    return max;
}
