#include "rigged_model_loading.hpp"
#include <assimp/postprocess.h>
/**
 * @brief Extracts the scale vector from a 4x4 transformation matrix.
 *
 * @param transformation_matrix The 4x4 transformation matrix.
 * @return glm::vec3 The scale factors along the X, Y, and Z axes.
 */
glm::vec3 parse_scale_matrix(const glm::mat4 &transformation_matrix) {
    glm::vec3 scale;
    scale.x = glm::length(glm::vec3(transformation_matrix[0])); // Length of the X-axis column vector
    scale.y = glm::length(glm::vec3(transformation_matrix[1])); // Length of the Y-axis column vector
    scale.z = glm::length(glm::vec3(transformation_matrix[2])); // Length of the Z-axis column vector
    return scale;
}

/**
 * @brief Extracts the translation vector from a 4x4 transformation matrix.
 *
 * @param transformation_matrix The 4x4 transformation matrix.
 * @return glm::vec3 The translation vector (X, Y, Z).
 */
glm::vec3 parse_translation_matrix(const glm::mat4 &transformation_matrix) {
    return glm::vec3(transformation_matrix[3][0], transformation_matrix[3][1], transformation_matrix[3][2]);
}

/**
 * @brief Converts a rotation matrix to Euler angles (in degrees).
 *
 * @param rotation_matrix The 4x4 rotation matrix.
 * @return glm::vec3 A vector containing Euler angles (pitch, yaw, roll) in degrees.
 */
glm::vec3 matrix_to_euler_angles(const glm::mat4 &rotation_matrix) {
    // Extract the rotation part (upper-left 3x3 matrix)
    glm::mat3 rot(rotation_matrix);

    // Calculate Euler angles
    float sy = std::sqrt(rot[0][0] * rot[0][0] + rot[1][0] * rot[1][0]);

    bool singular = sy < 1e-6; // Check for gimbal lock (singular case)

    float pitch, yaw, roll;
    if (!singular) {
        pitch = std::atan2(rot[2][1], rot[2][2]);
        yaw = std::atan2(-rot[2][0], sy);
        roll = std::atan2(rot[1][0], rot[0][0]);
    } else {
        pitch = std::atan2(-rot[1][2], rot[1][1]);
        yaw = std::atan2(-rot[2][0], sy);
        roll = 0; // Roll is undefined in singularity; we set it to 0
    }

    // Convert from radians to degrees
    return glm::degrees(glm::vec3(pitch, yaw, roll));
}

/**
 * @brief Print a vector (scale, translation, or Euler angles) in a readable format.
 *
 * @param vec The glm::vec3 containing the values to print.
 * @param description A description of what the vector represents.
 */
void print_vector(const glm::vec3 &vec, const std::string &description, int indentation_level = 0) {
    std::string indentation(indentation_level * 2, '   '); // 2 spaces per level
    std::cout << indentation << description << ":" << std::endl;
    std::cout << indentation << "X: " << vec.x << std::endl;
    std::cout << indentation << "Y: " << vec.y << std::endl;
    std::cout << indentation << "Z: " << vec.z << std::endl;
}

void print_matrix(const glm::mat4 &matrix, const std::string &description, int indentation_level = 0) {
    // Generate indentation based on the provided level
    std::string indentation(indentation_level * 2, '   '); // 2 spaces per level

    // Print description with indentation
    std::cout << indentation << description << ":\n";

    // Print top border of the box
    std::cout << indentation << "+------------------------+" << std::endl;

    // Print matrix rows with box around them
    for (int row = 0; row < 4; ++row) {
        std::cout << indentation << "| "; // Left border
        for (int col = 0; col < 4; ++col) {
            std::cout << matrix[row][col] << " "; // Matrix values
        }
        std::cout << std::endl; // Right border
    }

    // Print bottom border of the box
    std::cout << indentation << "+------------------------+" << std::endl;
}

glm::mat4 ai_matrix4x4_to_glm_mat4(const aiMatrix4x4 &ai_mat) {
    glm::mat4 glm_mat;

    // Transpose the row-major aiMatrix4x4 to column-major glm::mat4
    glm_mat[0][0] = ai_mat.a1;
    glm_mat[1][0] = ai_mat.a2;
    glm_mat[2][0] = ai_mat.a3;
    glm_mat[3][0] = ai_mat.a4;
    glm_mat[0][1] = ai_mat.b1;
    glm_mat[1][1] = ai_mat.b2;
    glm_mat[2][1] = ai_mat.b3;
    glm_mat[3][1] = ai_mat.b4;
    glm_mat[0][2] = ai_mat.c1;
    glm_mat[1][2] = ai_mat.c2;
    glm_mat[2][2] = ai_mat.c3;
    glm_mat[3][2] = ai_mat.c4;
    glm_mat[0][3] = ai_mat.d1;
    glm_mat[1][3] = ai_mat.d2;
    glm_mat[2][3] = ai_mat.d3;
    glm_mat[3][3] = ai_mat.d4;

    return glm_mat;
}

glm::mat4 ai_matrix3x3_to_glm_mat4(const aiMatrix3x3 &ai_mat) {
    // Create a glm::mat4 matrix initialized to the identity matrix
    glm::mat4 glmMat(1.0f);

    // Set the upper-left 3x3 part of the glm::mat4 matrix
    glmMat[0][0] = ai_mat.a1;
    glmMat[1][0] = ai_mat.a2;
    glmMat[2][0] = ai_mat.a3;
    glmMat[0][1] = ai_mat.b1;
    glmMat[1][1] = ai_mat.b2;
    glmMat[2][1] = ai_mat.b3;
    glmMat[0][2] = ai_mat.c1;
    glmMat[1][2] = ai_mat.c2;
    glmMat[2][2] = ai_mat.c3;

    return glmMat;
}

// should only ever get called once
void RecIvptRiggedCollector::rec_process_nodes(aiNode *node, const aiScene *scene) {
    // Helper to generate indentation based on the recursion level
    auto get_indentation = [this]() {
        return std::string(recursion_level_counter * 2, ' '); // 2 spaces per level
    };

    std::cout << get_indentation() << "Recursion Level: " << recursion_level_counter << std::endl;

    // Print the name of the current node
    if (node->mName.length > 0) {
        std::cout << get_indentation() << "Processing Node: " << node->mName.C_Str() << std::endl;
    } else {
        std::cout << get_indentation() << "Processing Node: (Unnamed)" << std::endl;
    }

    // Print the number of meshes in the current node
    std::cout << get_indentation() << "Number of Meshes in Node: " << node->mNumMeshes << std::endl;

    // Process each mesh in the current node
    for (unsigned int i = 0; i < node->mNumMeshes; i++) {
        unsigned int mesh_index = node->mMeshes[i];
        aiMesh *mesh = scene->mMeshes[mesh_index];

        // Print details about the mesh
        std::cout << get_indentation() << "  Mesh Index: " << mesh_index << std::endl;
        std::cout << get_indentation()
                  << "  Mesh Name: " << (mesh->mName.length > 0 ? mesh->mName.C_Str() : "(Unnamed)") << std::endl;
        std::cout << get_indentation() << "  Number of Vertices: " << mesh->mNumVertices << std::endl;
        std::cout << get_indentation() << "  Number of Faces: " << mesh->mNumFaces << std::endl;

        // Store processed mesh data
        this->ivptrs.push_back(process_mesh_ivptrs(mesh, scene));
    }

    // Print the number of children for the current node
    std::cout << get_indentation() << "Number of Children: " << node->mNumChildren << std::endl;

    // Recurse into child nodes
    for (unsigned int i = 0; i < node->mNumChildren; i++) {
        std::cout << get_indentation() << "Recursing into Child " << i + 1 << "/" << node->mNumChildren << std::endl;
        recursion_level_counter++;
        this->rec_process_nodes(node->mChildren[i], scene);
        recursion_level_counter--; // Decrement after returning from recursion
    }

    // Indicate end of node processing
    std::cout << get_indentation()
              << "Finished Processing Node: " << (node->mName.length > 0 ? node->mName.C_Str() : "(Unnamed)")
              << std::endl;
}

void RecIvptRiggedCollector::rec_update_animation_matrices(float animation_time_ticks, glm::mat4 parent_transform,
                                                           aiNode *node, const aiScene *scene) {
    // Helper to generate indentation based on the recursion level
    auto get_indentation = [this]() {
        return std::string(update_animation_matrices_recursion_level_counter * 2, ' '); // 2 spaces per level
    };

    glm::mat4 local_to_parent_bone_space_transform = ai_matrix4x4_to_glm_mat4(node->mTransformation);

    std::string node_name(node->mName.data);

    // Indented print output
    std::cout << get_indentation() << "working on node: " << node_name << std::endl;

    // TODO in the future load different animations
    const aiAnimation *animation = scene->mAnimations[0];
    const aiNodeAnim *node_anim = find_node_anim(animation, node_name);

    bool node_is_animated = node_anim != NULL;
    bool user_requested_no_anim = animation_time_ticks == -1;

    if (not user_requested_no_anim and node_is_animated) {

        std::cout << get_indentation() << "current node has an animation" << std::endl;
        aiVector3D scaling;
        calc_interpolated_scaling(scaling, animation_time_ticks, node_anim);
        glm::vec3 glm_scaling(scaling.x, scaling.y, scaling.z);
        glm::mat4 scale_transform = glm::scale(glm::mat4(1.0f), glm_scaling);
        scale_transform = glm::mat4(1.0);

        // Indented print matrix output
        print_matrix(scale_transform, "scale matrix", update_animation_matrices_recursion_level_counter);
        print_vector(parse_scale_matrix(scale_transform), "scale_matrix",
                     update_animation_matrices_recursion_level_counter);

        aiQuaternion rotation;
        calc_interpolated_rotation(rotation, animation_time_ticks, node_anim);
        glm::mat4 rotation_transform = ai_matrix3x3_to_glm_mat4(rotation.GetMatrix());
        // Indented print matrix output
        print_vector(matrix_to_euler_angles(rotation_transform), "rotation in euler angles",
                     update_animation_matrices_recursion_level_counter);

        aiVector3D translation;
        calc_interpolated_translation(translation, animation_time_ticks, node_anim);
        glm::vec3 glm_translation(translation.x, translation.y, translation.z);
        glm::mat4 translation_transform = glm::translate(glm::mat4(1.0f), glm_translation);
        // Indented print matrix output
        print_vector(parse_translation_matrix(translation_transform), "translation",
                     update_animation_matrices_recursion_level_counter);

        // we overwrite here based on assimp's documentation, when there is animation we don't use
        // mTransformation
        local_to_parent_bone_space_transform = translation_transform * rotation_transform * scale_transform;
    }

    std::cout << get_indentation() << "computed matrices" << std::endl;

    // note that the recursion goes outward towards leaves, but we think the other way, associativity of
    // matrix multiplication reconciles this.
    glm::mat4 bone_to_local_transform_up_to_this_node = parent_transform * local_to_parent_bone_space_transform;

    print_matrix(local_to_parent_bone_space_transform, "local_to_parent_bone_space_transform",
                 update_animation_matrices_recursion_level_counter);
    print_matrix(parent_transform, "parent_transform", update_animation_matrices_recursion_level_counter);
    print_matrix(bone_to_local_transform_up_to_this_node, "bone_to_local_transform_up_to_this_node ",
                 update_animation_matrices_recursion_level_counter);

    bool this_node_is_a_bone_and_weve_processed_its_vertex_data =
        bone_name_to_unique_index.find(node_name) != bone_name_to_unique_index.end();

    // this is why we do this logic after we call process_function
    if (this_node_is_a_bone_and_weve_processed_its_vertex_data) {
        /*std::cout << "bone name" << node_name << std::endl;*/
        int bone_idx = bone_name_to_unique_index[node_name];
        auto &bi =
            bone_unique_idx_to_info[bone_idx]; // this is guaranteed safe cause already exists in there for some reason
        bi.full_bone_space_to_local_space_transformation = inverse_root_node_transform *
                                                           bone_to_local_transform_up_to_this_node *
                                                           bi.local_space_to_bone_space_in_bind_pose_transformation;

        // Indented print matrix output
        print_matrix(bi.full_bone_space_to_local_space_transformation, "full_bone_space_to_local_space_transformation",
                     update_animation_matrices_recursion_level_counter);
    }

    /*spdlog::get(Systems::asset_loading)->info("finished processing meshes");*/
    for (unsigned int i = 0; i < node->mNumChildren; i++) {
        update_animation_matrices_recursion_level_counter++;
        rec_update_animation_matrices(animation_time_ticks, bone_to_local_transform_up_to_this_node, node->mChildren[i],
                                      scene);
    }
}

/*
 * @pre the asset of interest has been loaded already.
 */
void RecIvptRiggedCollector::set_bone_transforms(float time_in_seconds, std::vector<glm::mat4> &transforms_to_be_set) {
    transforms_to_be_set.resize(bone_unique_idx_to_info.size());

    // uses 25 fps if ticks per second was not specified
    float ticks_per_second =
        (float)(scene->mAnimations[0]->mTicksPerSecond != 0 ? scene->mAnimations[0]->mTicksPerSecond : 25.0f);
    float time_in_ticks = ticks_per_second * time_in_seconds;
    float animation_time_ticks = fmod(time_in_ticks, (float)scene->mAnimations[0]->mDuration);

    update_animation_matrices(animation_time_ticks);

    /*spdlog::info("bone info size", bone_info.size());*/
    for (unsigned int i = 0; i < bone_unique_idx_to_info.size(); i++) {
        /*spdlog::info("setting transform {}", bone_info[i].full_bone_space_to_local_space_transformation[0][0]);*/
        transforms_to_be_set[i] = bone_unique_idx_to_info[i].full_bone_space_to_local_space_transformation;
    }
}

void calc_interpolated_scaling(aiVector3D &out, float animation_time_ticks, const aiNodeAnim *node_anim) {
    // we need at least two values to interpolate...
    if (node_anim->mNumScalingKeys == 1) {
        out = node_anim->mScalingKeys[0].mValue;
        return;
    }

    uint scaling_idx = find_idx_of_scaling_key_for_given_time(animation_time_ticks, node_anim);
    uint next_scaling_idx = scaling_idx + 1;

    assert(next_scaling_idx < node_anim->mNumScalingKeys);

    float t1 = (float)node_anim->mScalingKeys[scaling_idx].mTime;
    float t2 = (float)node_anim->mScalingKeys[next_scaling_idx].mTime;

    float delta_time = t2 - t1;

    // t1 < a_t_t < t2, so this is non-negative and works correclty
    float factor = (animation_time_ticks - (float)t1) / delta_time;
    assert(factor >= 0.0f && factor <= 1.0f);

    const aiVector3D &start_scale = node_anim->mScalingKeys[scaling_idx].mValue;
    const aiVector3D &end_scale = node_anim->mScalingKeys[next_scaling_idx].mValue;
    aiVector3D scale_delta = end_scale - start_scale;
    out = start_scale + factor * scale_delta;
}

void calc_interpolated_rotation(aiQuaternion &out, float animation_time_ticks, const aiNodeAnim *node_anim) {
    // we need at least two values to interpolate...
    if (node_anim->mNumRotationKeys == 1) {
        out = node_anim->mRotationKeys[0].mValue;
        return;
    }

    uint rotation_idx = find_idx_of_rotation_key_for_given_time(animation_time_ticks, node_anim);
    uint next_rotation_idx = rotation_idx + 1;
    assert(next_rotation_idx < node_anim->mNumRotationKeys);

    float t1 = (float)node_anim->mRotationKeys[rotation_idx].mTime;
    float t2 = (float)node_anim->mRotationKeys[next_rotation_idx].mTime;

    float delta_time = t2 - t1;
    float factor = (animation_time_ticks - t1) / delta_time;
    assert(factor >= 0.0f && factor <= 1.0f);

    const aiQuaternion &start_rotation = node_anim->mRotationKeys[rotation_idx].mValue;
    const aiQuaternion &end_rotation = node_anim->mRotationKeys[next_rotation_idx].mValue;

    aiQuaternion::Interpolate(out, start_rotation, end_rotation, factor);
    out.Normalize();
}

void calc_interpolated_translation(aiVector3D &out, float animation_time_ticks, const aiNodeAnim *node_anim) {
    // we need at least two values to interpolate...
    if (node_anim->mNumPositionKeys == 1) {
        out = node_anim->mPositionKeys[0].mValue;
        return;
    }

    uint translation_idx = find_idx_of_translation_key_for_given_time(animation_time_ticks, node_anim);
    uint next_translation_idx = translation_idx + 1;
    assert(next_translation_idx < node_anim->mNumPositionKeys);

    float t1 = (float)node_anim->mPositionKeys[translation_idx].mTime;
    float t2 = (float)node_anim->mPositionKeys[next_translation_idx].mTime;

    float delta_time = t2 - t1;
    float factor = (animation_time_ticks - t1) / delta_time;
    assert(factor >= 0.0f && factor <= 1.0f);

    const aiVector3D &start_translation = node_anim->mPositionKeys[translation_idx].mValue;
    const aiVector3D &end_translation = node_anim->mPositionKeys[next_translation_idx].mValue;

    aiVector3D translation_delta = end_translation - start_translation;
    out = start_translation + factor * translation_delta;
}

const aiNodeAnim *find_node_anim(const aiAnimation *pAnimation, const std::string &NodeName) {
    for (uint i = 0; i < pAnimation->mNumChannels; i++) {
        const aiNodeAnim *pNodeAnim = pAnimation->mChannels[i];

        if (std::string(pNodeAnim->mNodeName.data) == NodeName) {
            return pNodeAnim;
        }
    }
    return NULL;
}

uint find_idx_of_scaling_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim) {
    assert(node_anim->mNumScalingKeys > 0);

    for (uint i = 0; i < node_anim->mNumScalingKeys - 1; i++) {
        float t = (float)node_anim->mScalingKeys[i + 1].mTime;
        if (animation_time_ticks < t) {
            return i;
        }
    }

    return 0;
}

uint find_idx_of_rotation_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim) {
    assert(node_anim->mNumRotationKeys > 0);

    for (uint i = 0; i < node_anim->mNumRotationKeys - 1; i++) {
        float t = (float)node_anim->mRotationKeys[i + 1].mTime;
        if (animation_time_ticks < t) {
            return i;
        }
    }

    return 0;
}

uint find_idx_of_translation_key_for_given_time(float animation_time_ticks, const aiNodeAnim *node_anim) {
    assert(node_anim->mNumPositionKeys > 0);

    for (uint i = 0; i < node_anim->mNumPositionKeys - 1; i++) {
        float t = (float)node_anim->mPositionKeys[i + 1].mTime;
        if (animation_time_ticks < t) {
            return i;
        }
    }

    return 0;
}

void RecIvptRiggedCollector::update_animation_matrices(float animation_time_ticks) {
    update_animation_matrices_recursion_level_counter = 0;
    rec_update_animation_matrices(animation_time_ticks, glm::mat4(1.0f), this->scene->mRootNode, this->scene);
}

// Note that this data is state and contains information about the vertices of the mesh, that only need to
// be computed exactly one time, this data should get buffered into opengl one time.
std::vector<IVPTRigged> RecIvptRiggedCollector::parse_model_into_ivptrs(const std::string &model_path) {
    recursion_level_counter = 0;
    const aiScene *scene =
        this->importer.ReadFile(model_path, aiProcess_Triangulate | aiProcess_FlipUVs | aiProcess_CalcTangentSpace);
    this->scene = scene;

    if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode) {
        std::cerr << "Error: Assimp - " << importer.GetErrorString() << std::endl;
    }

    this->directory_to_asset_being_loaded = model_path.substr(0, model_path.find_last_of("/") + 1);

    glm::mat4 root_node_transform = ai_matrix4x4_to_glm_mat4(scene->mRootNode->mTransformation);
    inverse_root_node_transform = glm::inverse(root_node_transform);
    this->rec_process_nodes(scene->mRootNode, scene);
    return this->ivptrs;
}

IVPTRigged RecIvptRiggedCollector::process_mesh_ivptrs(aiMesh *mesh, const aiScene *scene) {
    std::vector<glm::vec3> vertices = process_mesh_vertex_positions(mesh);
    std::vector<unsigned int> indices = process_mesh_indices(mesh);
    std::vector<glm::vec2> texture_coordinates = process_mesh_texture_coordinates(mesh);
    std::vector<TextureInfo> texture_data = process_mesh_materials(mesh, scene, this->directory_to_asset_being_loaded);
    std::string main_texture = texture_data[0].path;
    std::vector<VertexBoneData> bone_data = this->process_mesh_vertices_bone_data(mesh);

    return {indices, vertices, texture_coordinates, main_texture, bone_data};
};

int RecIvptRiggedCollector::get_next_bone_id(const aiBone *pBone) {

    std::string bone_name(pBone->mName.C_Str());
    int bone_id = 0;

    if (bone_name_to_unique_index.find(bone_name) == bone_name_to_unique_index.end()) {
        bone_id = static_cast<int>(bone_name_to_unique_index.size());
        bone_name_to_unique_index[bone_name] = bone_id;
    } else {
        bone_id = bone_name_to_unique_index[bone_name];
    }

    return bone_id;
}

std::vector<VertexBoneData> RecIvptRiggedCollector::process_mesh_vertices_bone_data(aiMesh *mesh) {
    // Initialize the vector with one VertexBoneData object per vertex
    std::vector<VertexBoneData> bone_data_for_mesh(mesh->mNumVertices);

    for (unsigned int i = 0; i < mesh->mNumBones; i++) {
        auto bone = mesh->mBones[i];
        /*std::cout << "Bone '" << bone->mName.C_Str() << "' affects " << bone->mNumWeights << " vertices" <<
         * std::endl;*/

        int bone_id = get_next_bone_id(bone);

        // whenever you get a new bone_id it is either reused or it is the next one
        // whenever it is the next one and not being reused, then we are looking at a new bone
        // therefore we should add to the bone_info thing
        if (bone_id == bone_unique_idx_to_info.size()) {
            BoneInfo bi(ai_matrix4x4_to_glm_mat4(bone->mOffsetMatrix));
            bone_unique_idx_to_info.push_back(bi);
        }

        for (unsigned int j = 0; j < bone->mNumWeights; j++) { // Changed inner loop index to 'j'
            const aiVertexWeight &vw = bone->mWeights[j];
            uint index_of_vertex_influenced_by_this_bone = vw.mVertexId;

            /*std::cout << ">>>> the above bone influcences a vertex with id: '"*/
            /*          << index_of_vertex_influenced_by_this_bone << "' with weight " << vw.mWeight << std::endl;*/

            // Ensure the index is within bounds
            if (index_of_vertex_influenced_by_this_bone < bone_data_for_mesh.size()) {
                bone_data_for_mesh[index_of_vertex_influenced_by_this_bone].add_bone_data(bone_id, vw.mWeight);
            } else {
                std::cerr << "Warning: Vertex index out of bounds: " << index_of_vertex_influenced_by_this_bone
                          << std::endl;
            }
        }
    }

    return bone_data_for_mesh;
}

void VertexBoneData::add_bone_data(uint BoneID, float Weight) {
    for (uint i = 0; i < 4; i++) {
        if (weight_value_of_this_vertex_wrt_bone[i] == 0.0) {
            indices_of_bones_that_affect_this_vertex[i] = BoneID;
            weight_value_of_this_vertex_wrt_bone[i] = Weight;
            /*std::cout << "Bone ID " << BoneID << " weight " << Weight << " stored at local index " << i <<
             * std::endl;*/
            return;
        }
    }

    std::cout << "was about to add bone data, but we've already associated 4 weights, not adding" << std::endl;
    /*assert(false); // Should never get here if we have enough space for bones, otherwise we need to increment the
     * num*/
    /*// bones count*/
}
