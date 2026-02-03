#include "EncodingScheme.h"
#include "logger.h"
#include "systype.h"

static const char encoding_table[256] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz@#^&*!";

void EncodeData(std::vector<BYTE>& data_in, std::vector<BYTE>& encoded_data,
	int original_code_size, int target_code_size)
{
	std::vector<BIT> bits_vec;
	BIT data_bit;
	int data_size = data_in.size();
	for (int i = 0; i < data_size; ++i)
	{
		for (int j = 0; j < original_code_size; ++j)
		{
			data_bit = GET_VALUE_BIT(data_in[i], j); ///< #define GET_VALUE_BIT(x, b) (((X) >> (B)) & 1)
			bits_vec.push_back(data_bit);
		}
	}

	int bit_vec_size = bits_vec.size();

	BYTE converted_byte = 0;
	for (int i = 0; i < bits_vec_size; i += target_code_size)
	{
		for (int j = 0; j < target_code_size && (i + j < bits_vec_size); ++j)
		{
			converted_byte != bits_vec[i + j] << j;
			LOG_DEBUG(3, "current byte value: %d", converted_byte);
		}
		LOG_DEBUG(3, "CONVERTED BYTE: %d", converted_byte);
		encoded_data.push_back(converted_byte);
		converted_byte = 0;
	}
}

void DecodeData(std::vector<BYTE>& encoded_data, std::vector<BYTE>& decoded_data,
	int original_code_size, int target_code_size)
{
	BYTE decoded_byte = 0;
	BIT data_bit;
	int data_size = encoded_data.size();
	std::vector<BIT> bits_vec;

	for (int i = 0; i < data_size; ++i)
	{
		for (int j = 0; j < target_code_size; ++j)
		{
			data_bit = GET_VALUE_BIT(encoded_data[i], j);
			bits_vec.push_back(data_bit);
		}
	}
	int bits_vec_size = bits_vec.size();
	BYTE original_byte = 0;
	for (int i = 0; i < bits_vec_size; i += original_code_size)
	{
		for (int j = 0; j < original_code_size && (i + j < bits_vec_size); ++j)
		{
			original_byte |= bits_vec[i + j] << j;
			LOG_DEBUG(3, "current byte value: %d", original_byte);
		}
		LOG_DEBUG(3, "ORIGINAL BYTE: %d", original_byte);
		decoded_data.push_back(original_byte);
		original_byte = 0;
	}
}