/* Copyright (C) 2026 JLinAlg contributors
 * SPDX-License-Identifier: GPL-2.0-or-later */
#include <jni.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <cholmod.h>

typedef struct
{
    cholmod_common common;
    cholmod_sparse *matrix;
    cholmod_factor *factor;
    int dimension;
} jlinalg_cholmod_factor;

static void throw_java(JNIEnv *env, const char *type, const char *message)
{
    jclass exception = (*env)->FindClass(env, type);
    if (exception != NULL) (*env)->ThrowNew(env, exception, message);
}

static jlinalg_cholmod_factor *from_handle(jlong value)
{
    return (jlinalg_cholmod_factor *) (uintptr_t) value;
}

static int numeric_factorize(jlinalg_cholmod_factor *handle)
{
    handle->common.status = CHOLMOD_OK;
    return cholmod_factorize(handle->matrix, handle->factor, &handle->common)
        && handle->common.status >= CHOLMOD_OK
        && handle->factor->minor == (size_t) handle->dimension;
}

static void destroy_factor(jlinalg_cholmod_factor *handle)
{
    if (handle == NULL) return;
    if (handle->factor != NULL)
        cholmod_free_factor(&handle->factor, &handle->common);
    if (handle->matrix != NULL)
        cholmod_free_sparse(&handle->matrix, &handle->common);
    cholmod_finish(&handle->common);
    free(handle);
}

JNIEXPORT jlong JNICALL
Java_org_jlinalg_compute_CholmodNative_create(
    JNIEnv *env, jclass type, jint dimension, jintArray row_starts,
    jintArray column_indices, jdoubleArray values, jboolean lower,
    jboolean natural_ordering)
{
    (void) type;
    jsize nonzeros = (*env)->GetArrayLength(env, values);
    jlinalg_cholmod_factor *handle = calloc(1, sizeof(*handle));
    if (handle == NULL)
    {
        throw_java(env, "java/lang/OutOfMemoryError", "cannot allocate CHOLMOD handle");
        return 0;
    }
    if (!cholmod_start(&handle->common))
    {
        free(handle);
        throw_java(env, "java/lang/IllegalStateException", "cannot initialize CHOLMOD");
        return 0;
    }
    handle->dimension = dimension;
    handle->common.final_ll = 1;
    handle->common.supernodal = CHOLMOD_AUTO;
    handle->common.nmethods = 1;
    handle->common.method[0].ordering = natural_ordering
        ? CHOLMOD_NATURAL : CHOLMOD_AMD;
    handle->common.postorder = natural_ordering ? 0 : 1;
    handle->matrix = cholmod_allocate_sparse(
        (size_t) dimension, (size_t) dimension, (size_t) nonzeros,
        1, 1, lower ? 1 : -1, CHOLMOD_REAL, &handle->common);
    if (handle->matrix == NULL)
    {
        destroy_factor(handle);
        throw_java(env, "java/lang/OutOfMemoryError", "cannot allocate CHOLMOD matrix");
        return 0;
    }
    (*env)->GetIntArrayRegion(env, row_starts, 0, dimension + 1,
        (jint *) handle->matrix->p);
    (*env)->GetIntArrayRegion(env, column_indices, 0, nonzeros,
        (jint *) handle->matrix->i);
    (*env)->GetDoubleArrayRegion(env, values, 0, nonzeros,
        (jdouble *) handle->matrix->x);
    if ((*env)->ExceptionCheck(env))
    {
        destroy_factor(handle);
        return 0;
    }
    /* JDistlib CSR uses one-based MKL indexing. CHOLMOD uses zero-based
       compressed-column indexing; interpreting CSR(A) as CSC(A') also lets
       us pass the stored triangle without an explicit transpose. */
    int32_t *column_starts = handle->matrix->p;
    int32_t *row_indices = handle->matrix->i;
    for (int column = 0; column <= dimension; column++)
        column_starts[column]--;
    for (jsize index = 0; index < nonzeros; index++)
        row_indices[index]--;
    handle->factor = cholmod_analyze(handle->matrix, &handle->common);
    if (handle->factor == NULL || !numeric_factorize(handle))
    {
        destroy_factor(handle);
        throw_java(env, "java/lang/IllegalArgumentException",
            "CHOLMOD matrix is not symmetric positive definite");
        return 0;
    }
    return (jlong) (uintptr_t) handle;
}

JNIEXPORT void JNICALL
Java_org_jlinalg_compute_CholmodNative_refactor(
    JNIEnv *env, jclass type, jlong native_handle, jdoubleArray values)
{
    (void) type;
    jlinalg_cholmod_factor *handle = from_handle(native_handle);
    if (handle == NULL)
    {
        throw_java(env, "java/lang/IllegalStateException", "CHOLMOD factor is closed");
        return;
    }
    jsize nonzeros = (*env)->GetArrayLength(env, values);
    if ((size_t) nonzeros != handle->matrix->nzmax)
    {
        throw_java(env, "java/lang/IllegalArgumentException",
            "CHOLMOD refactor values have the wrong length");
        return;
    }
    (*env)->GetDoubleArrayRegion(env, values, 0, nonzeros,
        (jdouble *) handle->matrix->x);
    if ((*env)->ExceptionCheck(env)) return;
    if (!numeric_factorize(handle))
        throw_java(env, "java/lang/IllegalArgumentException",
            "CHOLMOD matrix is not symmetric positive definite");
}

JNIEXPORT jdouble JNICALL
Java_org_jlinalg_compute_CholmodNative_logDeterminant(
    JNIEnv *env, jclass type, jlong native_handle)
{
    (void) type;
    jlinalg_cholmod_factor *handle = from_handle(native_handle);
    if (handle == NULL)
    {
        throw_java(env, "java/lang/IllegalStateException", "CHOLMOD factor is closed");
        return NAN;
    }
    cholmod_factor *factor = handle->factor;
    double result = 0.0;
    if (factor->is_super)
    {
        int32_t *super = factor->super;
        int32_t *pi = factor->pi;
        int32_t *px = factor->px;
        double *x = factor->x;
        for (size_t node = 0; node < factor->nsuper; node++)
        {
            int32_t columns = super[node + 1] - super[node];
            int32_t leading = pi[node + 1] - pi[node];
            for (int32_t column = 0; column < columns; column++)
            {
                double diagonal = x[px[node] + column * leading + column];
                if (!(diagonal > 0.0))
                {
                    throw_java(env, "java/lang/IllegalStateException",
                        "CHOLMOD factor has a non-positive diagonal");
                    return NAN;
                }
                result += 2.0 * log(diagonal);
            }
        }
    }
    else
    {
        int32_t *p = factor->p;
        double *x = factor->x;
        for (int column = 0; column < handle->dimension; column++)
        {
            double diagonal = x[p[column]];
            if (!(diagonal > 0.0))
            {
                throw_java(env, "java/lang/IllegalStateException",
                    "CHOLMOD factor has a non-positive diagonal");
                return NAN;
            }
            result += (factor->is_ll ? 2.0 : 1.0) * log(diagonal);
        }
    }
    return result;
}

JNIEXPORT jint JNICALL
Java_org_jlinalg_compute_CholmodNative_factorNonzeroCount(
    JNIEnv *env, jclass type, jlong native_handle)
{
    (void) type;
    jlinalg_cholmod_factor *handle = from_handle(native_handle);
    if (handle == NULL)
    {
        throw_java(env, "java/lang/IllegalStateException", "CHOLMOD factor is closed");
        return 0;
    }
    int32_t *counts = handle->factor->ColCount;
    int64_t result = 0;
    for (int column = 0; column < handle->dimension; column++)
        result += counts[column];
    return result > INT32_MAX ? INT32_MAX : (jint) result;
}

JNIEXPORT jintArray JNICALL
Java_org_jlinalg_compute_CholmodNative_permutation(
    JNIEnv *env, jclass type, jlong native_handle)
{
    (void) type;
    jlinalg_cholmod_factor *handle = from_handle(native_handle);
    if (handle == NULL)
    {
        throw_java(env, "java/lang/IllegalStateException", "CHOLMOD factor is closed");
        return NULL;
    }
    jintArray result = (*env)->NewIntArray(env, handle->dimension);
    if (result != NULL)
        (*env)->SetIntArrayRegion(env, result, 0, handle->dimension,
            (jint *) handle->factor->Perm);
    return result;
}

JNIEXPORT void JNICALL
Java_org_jlinalg_compute_CholmodNative_solveInPlace(
    JNIEnv *env, jclass type, jlong native_handle, jdoubleArray right_hand_side,
    jint right_hand_sides)
{
    (void) type;
    jlinalg_cholmod_factor *handle = from_handle(native_handle);
    if (handle == NULL)
    {
        throw_java(env, "java/lang/IllegalStateException", "CHOLMOD factor is closed");
        return;
    }
    cholmod_dense *input = cholmod_allocate_dense(
        handle->dimension, right_hand_sides, handle->dimension,
        CHOLMOD_REAL, &handle->common);
    if (input == NULL)
    {
        throw_java(env, "java/lang/OutOfMemoryError", "cannot allocate CHOLMOD right-hand side");
        return;
    }
    jdouble *java_values = (*env)->GetDoubleArrayElements(
        env, right_hand_side, NULL);
    if (java_values == NULL)
    {
        cholmod_free_dense(&input, &handle->common);
        return;
    }
    double *native_values = input->x;
    for (int row = 0; row < handle->dimension; row++)
        for (int column = 0; column < right_hand_sides; column++)
            native_values[column * handle->dimension + row] =
                java_values[row * right_hand_sides + column];
    cholmod_dense *solution = cholmod_solve(
        CHOLMOD_A, handle->factor, input, &handle->common);
    int solved_ok = solution != NULL;
    if (solution != NULL)
    {
        double *solved = solution->x;
        for (int row = 0; row < handle->dimension; row++)
            for (int column = 0; column < right_hand_sides; column++)
                java_values[row * right_hand_sides + column] =
                    solved[column * handle->dimension + row];
    }
    (*env)->ReleaseDoubleArrayElements(env, right_hand_side, java_values, 0);
    cholmod_free_dense(&solution, &handle->common);
    cholmod_free_dense(&input, &handle->common);
    if (!solved_ok)
        throw_java(env, "java/lang/IllegalStateException", "CHOLMOD solve failed");
}

JNIEXPORT void JNICALL
Java_org_jlinalg_compute_CholmodNative_destroy(
    JNIEnv *env, jclass type, jlong native_handle)
{
    (void) env;
    (void) type;
    destroy_factor(from_handle(native_handle));
}

JNIEXPORT jstring JNICALL
Java_org_jlinalg_compute_CholmodNative_version(JNIEnv *env, jclass type)
{
    (void) type;
    return (*env)->NewStringUTF(env, JLINALG_CHOLMOD_VERSION);
}
