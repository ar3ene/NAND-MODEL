#pragma once

#ifndef _ENCODING_SCHEME_H_
#define _ENCODING_SCHEME_H_

#include <Windows.h>
#include <vector>

void EncodeData(std::vector<BYTE>& data_in, std::vector<BYTE>& encoded_data,
	int original_code_size, int target_code_size);
void DecodeData(std::vector<BYTE>& encoded_data, std::vector<BYTE>& decoded_data,
	int original_code_size, int target_code_size);

#endif // ! _ENCODING_SCHEME_H_